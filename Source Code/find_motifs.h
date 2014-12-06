#ifndef _FIND_MOTIFS_H_
#define _FIND_MOTIFS_H_

#include "common.h"
#include "candidate.h"
#include "subsequence.h"
#include "ucr_dtw.h"
#include "trivial_match_map.h"
#include <limits>
#include <fstream>
#include <queue>


bool operator<(const Candidate& lhs, const Candidate& rhs) {
    //this allows us to use candidates as the key in std::map for looking up winner tree external indices
    return lhs.query_loc < rhs.query_loc || (lhs.query_loc == rhs.query_loc && lhs.loc < rhs.loc);
}

bool operator>(const Candidate& lhs, const Candidate& rhs) {
    //this is used for outputting top candidates in order
    return lhs.query_loc > rhs.query_loc || (lhs.query_loc == rhs.query_loc && lhs.loc > rhs.loc);
}

//the following class uses a max winner tree to hold the top-K items, with weaker matches having a higher priority
//and consequentially allowing weakest match lookup in constant time and replacing the weakest match or an arbitrary
//candidate in logarithmic time
class TopCandidates {
private:
    struct Node {
        Candidate candidate;
        //is_placeholder indicates whether this holds a real candidate or is just an array slot placeholder
        bool is_placeholder;
        bool contains_candidate;
        void reset() {
            candidate.reset();
            contains_candidate = false;
        }
        bool operator==(const Node& rhs) const {
            return is_placeholder == rhs.is_placeholder || candidate == rhs.candidate;
        }
    };
    std::map<Candidate, size_t> winner_tree_external_indices;

    //note: packed array is 1-based to make the math easier
    std::vector<Node> tree;
    size_t num_nodes;
    size_t K;
    size_t num_stored_candidates = 0;

    size_t get_external_index_by_candidate(Candidate const& candidate) {
        if (winner_tree_external_indices.find(candidate) == winner_tree_external_indices.end())
            //candidate not found
            //xxx consider not throwing an exception here, since that could interupt a long-running process
            throw std::invalid_argument("Candidate not found in external index lookup map");
        else
            return winner_tree_external_indices[candidate];
    }
    /*
        start debug stuff
    */
        size_t num_children(size_t subtree_root_index) const {
            size_t n = 0;
            if (subtree_root_index*2 > num_nodes) return n;
            ++n;
            if (subtree_root_index*2 + 1 > num_nodes) return n;
            ++n;
            n += num_children(subtree_root_index*2);
            n += num_children(subtree_root_index*2 + 1);
            return n;
        }
        void write_subtree_buffer(size_t subtree_root_index,
                                  std::string buffer_lines[],
                                  size_t root_line_index,
                                  size_t lbound_line_index /*inclusive*/,
                                  size_t ubound_line_index /*exclusive*/)
        {
            Node subtree_root = tree[subtree_root_index];
            std::ostringstream oss;
            //print the node
            Candidate c = subtree_root.candidate;
            if (subtree_root.is_placeholder)
                oss << "\x1b[31m";
            size_t external_index = 0;
            try {
                external_index = get_external_index_by_candidate(c);
            } catch(...) {}
            oss << "[" << subtree_root_index << "]: [Q: " << c.query_loc << ", C: " << c.loc << ", D: " << c.dist << ", EI: " << external_index << "]";
            if (subtree_root.is_placeholder)
                oss << "\x1b[0m";
            buffer_lines[root_line_index] += oss.str();
            //print the right descendents
            size_t right_index = subtree_root_index*2 + 1;
            if (right_index <= num_nodes) {
                //at least 1 right child
                size_t top_dashes = 1;
                Node const& right_child = tree[right_index];
                size_t left_child_index = right_index * 2;
                if (left_child_index <= num_nodes) {
                    //right child has at least 1 left child
                    top_dashes += 2 * (1 + num_children(left_child_index));
                }
                size_t top_line_index = root_line_index - 1;
                while (top_line_index >= root_line_index - top_dashes)
                    buffer_lines[top_line_index--] += "|  ";
                size_t right_child_line_index = top_line_index;
                buffer_lines[top_line_index--] += "+--";
                while (top_line_index >= lbound_line_index)
                    buffer_lines[top_line_index--] += "   ";
                write_subtree_buffer(right_index,
                                     buffer_lines,
                                     right_child_line_index,
                                     lbound_line_index,
                                     root_line_index);
            }
            //print the left descendents
            size_t left_index = subtree_root_index*2;
            if (left_index <= num_nodes) {
                //at least 1 left child
                size_t bottom_dashes = 1;
                Node const& left_child = tree[left_index];
                size_t right_child_index = right_index * 2 + 1;
                if (right_child_index <= num_nodes) {
                    //right child has at least 1 left child
                    bottom_dashes += 2 * (1 + num_children(right_child_index));
                }
                size_t bottom_line_index = root_line_index + 1;
                while (bottom_line_index <= root_line_index + bottom_dashes)
                    buffer_lines[bottom_line_index++] += "|  ";
                size_t left_child_line_index = bottom_line_index;
                buffer_lines[bottom_line_index++] += "+--";
                while (bottom_line_index < ubound_line_index)
                    buffer_lines[bottom_line_index++] += "   ";
                write_subtree_buffer(left_index,
                                     buffer_lines,
                                     left_child_line_index,
                                     root_line_index + 1,
                                     ubound_line_index);
            }
        }
        std::ostream& print(std::ostream& out) {
            size_t root_index = 1;
            size_t num_lines = num_nodes * 2 - 1;
            std::string buffer_lines[num_lines + 1];
            size_t root_line_index = 1;
            size_t right_index = 1 * 2 + 1;
            if (right_index <= num_nodes) {
                root_line_index += 2 * (1 + num_children(right_index));
            }
            write_subtree_buffer(root_index, buffer_lines, root_line_index, 1, num_lines + 1);
            for (size_t i = 1; i <= num_lines; ++i)
                out << buffer_lines[i] << std::endl;
            return out;
        }

    /*
        end debug stuff
    */


    size_t prepare_subtree(size_t index) {
        //since all slots represent candidates which are initially the same distance, we just need to
        //bubble up stored external indices (shaves lg(N) comparisons off replacing the weakest candidate
        //and O(N) off of eliminating an arbitrary candidate) and whether the node is a placeholder. this
        //function is recursive and returns an index to an external node slot. we traverse down the left child, so
        //this function returns the index to the external node arrived at if we were to follow the left leg
        //from the current position to the bottom of the tree. consequentially, adding items preserves tree
        //completeness, and we don't need manually ensure that we're storing no more than K items
        size_t external_index;
        size_t left_child_index = index * 2;
        size_t right_child_index = left_child_index + 1;
        bool does_left_child_exist = (left_child_index <= num_nodes);
        if (does_left_child_exist) {
            external_index = prepare_subtree(left_child_index);
            //right child guaranteed to exist if left child exists
            prepare_subtree(right_child_index);
        }
        else {
            //this is an external node, so return the current index, which may or may not be relevant to the caller
            external_index = index;
        }
        //the following has the effect of bubbling placeholder designations up to where the left sibling is a non-placeholder
        tree[index] = tree[external_index];
        return external_index;
    }
    Node const& get_winning_child(size_t index) const {
        //xxx this does *not* do bounds checking. if index points to an external node this will exhibit undefined behavior
        const Node& left_child = tree[2 * index];
        const Node& right_child = tree[2 * index + 1];
        if (right_child.is_placeholder == false && right_child.candidate.dist > left_child.candidate.dist)
            return right_child;
        else
            return left_child;
    }
    size_t recursive_replace_weakest_match(Candidate const& new_candidate, size_t subtree_root_index) {
        //replaces the relevant external node with the new candidate and replays up the tree, returning the
        //index of the replaced external node
        Node& subtree_root = tree[subtree_root_index];
        size_t left_child_index = 2 * subtree_root_index;
        size_t right_child_index = left_child_index + 1;
        size_t external_index;
        if (left_child_index <= num_nodes) {
            //subtree root is internal
            if (subtree_root.candidate == tree[left_child_index].candidate)
                external_index = recursive_replace_weakest_match(new_candidate, left_child_index);
            else
                //xxx: we assume the right child will never be a placeholder - is this valid?
                external_index = recursive_replace_weakest_match(new_candidate, right_child_index);
            subtree_root = get_winning_child(subtree_root_index);
        } else {
            //subtree root is external
            subtree_root.candidate = new_candidate;
            external_index = subtree_root_index;
            if ( ! subtree_root.contains_candidate) {
                //this is the only time we might be adding a candidate to the tree
                subtree_root.contains_candidate = true;
                ++num_stored_candidates;
            }
        }
        return external_index;
    }
    void replay_up(size_t external_index) {
        for(size_t i = external_index / 2; i != 0; i /= 2) {
            Node winner = get_winning_child(i);
            if (tree[i] == winner)
                break; //contest outcome didn't change
            tree[i] = winner;
        }
    }
    size_t count_internal_nodes() const {
        return (num_nodes - 1) / 2;
    }
public:
    class TopCandidates_Iter: public std::iterator<std::forward_iterator_tag, Candidate> {
    private:
        std::vector<Node>::iterator tree_iter;
        typedef TopCandidates_Iter self_type;
    public:
        TopCandidates_Iter(self_type const& src): tree_iter(src.tree_iter) {}
        explicit TopCandidates_Iter(std::vector<Node>* tree, size_t start_index): tree_iter(tree->begin() + start_index) {}
        reference operator*() const {
            return tree_iter->candidate;
        }
        pointer operator->() const {
            return &operator*();
        }
        self_type& operator=(self_type const& src) = delete;
        self_type& operator++() {
            //preincrement
            ++tree_iter;
            return *this;
        }
        self_type operator++(int) = delete;
        bool operator==(TopCandidates_Iter const& rhs) const {
            return tree_iter == rhs.tree_iter;
        }
        bool operator!=(TopCandidates_Iter const& rhs) const {
            return ! operator==(rhs);
        }
    };
    typedef TopCandidates_Iter iterator;
    iterator begin() {
        size_t num_internal_nodes = count_internal_nodes();
        size_t first_external_index = num_internal_nodes + 1;
        iterator iter(&tree, first_external_index);
        return iter;
    }
    iterator end() {
        size_t num_internal_nodes = count_internal_nodes();
        size_t first_external_index = num_internal_nodes + 1;
        size_t first_unfilled_node = first_external_index + num_stored_candidates;
        iterator iter(&tree, first_unfilled_node);
        return iter;
    }
    TopCandidates(size_t max_entries): K(max_entries) {
        //we need 2^n external nodes for easy math, where n is an integer
        //this way, if we're at an internal node, we can always assume the right child exists
        size_t num_external_nodes = pow(2, ceil(lg(K)));
        size_t num_internal_nodes = num_external_nodes - 1;
        num_nodes = num_external_nodes + num_internal_nodes;
        tree.resize(num_nodes + 1);
        tree[0].is_placeholder = true;
        //since we can only store K candidates in the tree, disable external nodes more than K slots from the left.
        //doing so automagically guarantees we store at most K items in the tree
        for (unsigned int i = num_nodes; i > num_internal_nodes + K; --i)
            tree[i].is_placeholder = true;
        prepare_subtree(1);
        //print(std::cout) << std::endl;
    }
    Candidate const& weakest_candidate() const {
        return tree[1].candidate;
    }
    size_t weakest_candidate_external_index() {
        const Candidate& c = weakest_candidate();
        return get_external_index_by_candidate(c);
    }
    void replace_candidate(Candidate const& old_candidate, Candidate const& new_candidate) {
        //std::cout << "before:" << std::endl;
        //print(std::cout) << std::endl;
        size_t external_index = get_external_index_by_candidate(old_candidate);
        tree[external_index].candidate = new_candidate;
        winner_tree_external_indices.erase(old_candidate);
        winner_tree_external_indices[new_candidate] = external_index;
        replay_up(external_index);
        //std::cout << "after:" << std::endl;
        //print(std::cout) << std::endl;
    }
    void replace_weakest_candidate(Candidate const& new_candidate) {
        //std::cout << "before:" << std::endl;
        //print(std::cout) << std::endl;
        size_t external_index = recursive_replace_weakest_match(new_candidate, 1);
        winner_tree_external_indices[new_candidate] = external_index;
        //std::cout << "after:" << std::endl;
        //print(std::cout) << std::endl;
    }
    void reset() {
        //this function allows us to reuse this class instance for a new query.
        //reset the left K external nodes
        iterator iter_end = end();
        for (iterator iter = begin(); iter != iter_end; ++iter) {
            size_t external_index = get_external_index_by_candidate(*iter);
            tree[external_index].reset();
        }
        //bubble the reset nodes up the tree
        prepare_subtree(1);
        //print(std::cout) << std::endl;
        num_stored_candidates = 0;
    }
};

class MotifFinder {
private:
    std::vector<double> time_series;
    TrivialMatchMap trivial_matches_map;
    std::ofstream results_output;
    const char field_separator;
    TopCandidates top_k;
    SubsequenceLookup subsequences;
    size_t curr_query_loc;
    void output_and_disregard_top_k() {
        TopCandidates::iterator iter_end = top_k.end();
        //use a priority queue so that we output the candidates in order of candidate location (important for post-processing)
        std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> output_pq(top_k.begin(), top_k.end());
        while ( ! output_pq.empty()) {
            Candidate const& c = output_pq.top();
            results_output << c.query_loc << field_separator << c.loc << field_separator << c.dist << std::endl;
            output_pq.pop();
        }
        top_k.reset();
        trivial_matches_map.reset();
    }
public:
    MotifFinder(std::string const& series_filepath, size_t series_len, unsigned int K, std::string const& output_filepath, const char results_field_separator):
        time_series(series_len),
        results_output(output_filepath, std::ofstream::out | std::ofstream::trunc),
        field_separator(results_field_separator),
        top_k(K),
        subsequences(time_series)
    {
        std::cerr << "Reading time series data file" << std::endl;
        std::ifstream in(series_filepath);
        size_t i = 0;
        while (in >> time_series[i++]);
        std::cerr << "first entry: " << time_series[0] << std::endl;
        in.close();
    }
    double weakest_distance() {
        const Candidate& weakest_candidate = top_k.weakest_candidate();
        return weakest_candidate.dist;
    }
    void try_register_candidate(size_t new_candidate_index, double distance) {
        //note: this function does not check to see that the new candidate is actually better and worth remembering.
        //that responsibility is left to the client, who should be checking that with weakest_distance
        Candidate new_candidate;
        new_candidate.dist = distance;
        new_candidate.query_loc = curr_query_loc;
        new_candidate.loc = new_candidate_index;
        Candidate old_candidate;
        if (trivial_matches_map.detect_trivial_match(new_candidate, &old_candidate)) {
            //found a trivial match
            //the following assert assumes we want to reset after incrementing the query
            assert(old_candidate.query_loc == new_candidate.query_loc);
            if (new_candidate.dist < old_candidate.dist) {
                //old candidate is trivial match (new is better), so disregard old and register the new one
                trivial_matches_map.push_candidate(new_candidate);
                top_k.replace_candidate(old_candidate, new_candidate);
            }
        } else {
            //no trivial match
            trivial_matches_map.push_candidate(new_candidate);
            top_k.replace_weakest_candidate(new_candidate);
        }
    }
    void run(size_t candidate_start_offset, size_t candidate_increment) {
        using namespace std::placeholders;
        UCR_DTW subsequence_search(subsequences);
        auto weakest_dist_callback = std::bind(&MotifFinder::weakest_distance, this);
        auto register_candidate_callback = std::bind(&MotifFinder::try_register_candidate, this, _1, _2);
        for (size_t i = 0; i < TIME_SERIES_LEN - QUERY_LEN; ++i) {
            curr_query_loc = i;
            subsequence_search.single_pass(i, candidate_start_offset, candidate_increment, weakest_dist_callback, register_candidate_callback);
            output_and_disregard_top_k();
        }
    }
    std::vector<double> const& get_timeseries() {
        return time_series;
    }
    SubsequenceLookup& get_subsequence_lookup() {
        return subsequences;
    }
};

#endif
