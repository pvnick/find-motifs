#ifndef _FIND_MOTIFS_H_
#define _FIND_MOTIFS_H_

#include "common.h"
#include "cache.h"
#include "ucr_dtw.h"
#include <limits>

class TopCandidates;

struct Candidate {
friend class TopCandidates;
public:
    typedef double DIST_TYPE;
private:
    static constexpr double max_dist = std::numeric_limits<DIST_TYPE>::max();
    size_t winner_tree_external_index;
public:
    DIST_TYPE dist;
    size_t loc;
    size_t query_loc;
    //if the candidate held in this node is, at a subsequent query position, determined to be a trivial match,
    //is_eliminated will get set to true, and this candidate will be skipped when candidates are outputted
    bool is_eliminated;
    Candidate() {
        reset();
    };
    void eliminate() {
        reset();
        is_eliminated = true;
    }
    bool is_valid() const {
        //since candidates must be at least QUERY_LEN points past the start of the time series, their location can't be zero
        return loc != 0;
    }
    void reset() {
        is_eliminated = false;
        dist = max_dist;
        loc = 0;
        query_loc = 0;
    }
};

//the following class uses a max winner tree to hold the top-K items, with weaker matches having a higher priority
//and consequentially allowing weakest match lookup in constant time and replacing the weakest match or an arbitrary
//candidate in logarithmic time
class TopCandidates {
private:
    struct Node {
        Candidate candidate;
        //is_placeholder indicates whether this holds a real candidate or is just an array slot placeholder
        bool is_placeholder;
        void reset() {
            candidate.reset();
        }
    };
    //note: packed array is 1-based to make the math easier
    Node* tree;
    size_t num_nodes;
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
        bool is_left_internal = (left_child_index <= num_nodes / 2);
        if (is_left_internal) {
            //if left child is an internal node, then right child is guaranteed to be internal
            external_index = prepare_subtree(left_child_index);
            size_t right_child_index = left_child_index + 1;
            prepare_subtree(right_child_index);
        } else {
            //this is an external node, so return the current index, which may or may not be relevant to the caller
            external_index = index;
        }
        //the following has the effect of bubbling placeholder designations up to where the left sibling is a non-placeholder
        tree[index] = tree[external_index]; 
        return external_index;
    }
    Node get_winning_child(size_t index) const {
        //xxx this does *not* do bounds checking. if index points to an external node this will display undefined behavior
        Node left_child = tree[2 * index];
        Node right_child = tree[2 * index + 1];
        if (right_child.is_placeholder == false && right_child.candidate.dist > left_child.candidate.dist)
            return right_child;
        else
            return left_child;
    }
    void replay_up(size_t external_index) {
        for(size_t i = external_index / 2; i != 0; i /= 2) {
            Node winner = get_winning_child(i);
            if (tree[i] == winner)
                break; //contest outcome didn't change
            tree[i] = winner;
        }
    }
    size_t count_internal_nodes() {
        return (num_nodes - 1) / 2;
    }
public:
    class TopCandidates_Iter: public std::iterator<std::forward_iterator_tag, Candidate> {
    //we don't care about order when iterating through the top-k matches, so just traverse the winner tree external nodes
    friend class TopCandidates;
    private:
        typedef TopCandidates_Iter self_type;
        value_type* tree;
        size_t curr_index;
        //end_index is typically set to one past the last occupiable tree node (ie the first placeholder 
        //external node in a non-full tree)
        size_t end_index;
        //no reason to expose the normal constructor to the client of TopCandidates
        explicit TopCandidates_Iter(value_type* tree_array, size_t tree_array_index, size_t max_index): 
            tree(tree_array), 
            curr_index(tree_array_index),
            end_index(max_index) {}
    public:
        TopCandidates_Iter(const self_type& src): tree(src.tree), curr_index(src.curr_index) {}
        reference operator*() const {
            return tree[curr_index].candidate;
        }
        pointer operator->() const {
            return &operator*();
        }
        self_type& operator=(const self_type& src) = delete;
        self_type& operator++() {
            //preincrement
            ++curr_index;
            return *this;
        }
        self_type operator++(int) = delete;
        bool operator==(const TopCandidates_Iter& rhs) const {
            return curr_index == rhs.curr_index;
        }
        bool operator!=(const TopCandidates_Iter& rhs) const {
            return ! operator==(rhs);
        }
    };
    typedef TopCandidates_Iter iterator;
    iterator begin() {
        size_t num_internal_nodes = count_internal_nodes();
        size_t first_external_index = num_internal_nodes + 1;
        size_t first_external_placeholder = first_external_index + K;
        iterator iter(tree, first_external_index, first_external_placeholder);
        return iter;
    }
    iterator end() {
        size_t num_internal_nodes = count_internal_nodes();
        size_t first_external_index = num_internal_nodes + 1;
        size_t first_external_placeholder = first_external_index + K;
        iterator iter(tree, first_external_placeholder, first_external_placeholder);
        return iter;
    }
    TopCandidates(unsigned int K) {
        //we need 2^n external nodes for easy math, where n is an integer
        size_t num_external_nodes = pow(2, ceil(lg(K)));
        size_t num_internal_nodes = num_external_nodes - 1;
        num_nodes = num_external_nodes + num_internal_nodes;
        tree = new Node[num_nodes];
        tree[0].is_placeholder = true;
        //since we can only store K candidates in the tree, disable external nodes more than K slots from the left.
        //doing so automagically guarantees we store at most K items in the tree
        for (int i = num_nodes; i > num_internal_nodes + K; --i)
            tree[i].is_placeholder = true;
        prepare_subtree(1);
    }
    ~TopCandidates() {
        delete[] tree;
    }
    Candidate weakest_candidate() const {
        return tree[1].candidate;
    }
    size_t weakest_candidate_external_index() const {
        Candidate m = weakest_candidate();
        return m.winner_tree_external_index;
    }
    void replace_candidate(Candidate old_candidate, Candidate new_candidate) {
        size_t external_index = old_candidate.winner_tree_external_index;
        new_candidate.winner_tree_external_index = external_index;
        tree[external_index].candidate = new_candidate;
        replay_up(external_index);
    }
    void replace_weakest_candidate(Candidate new_candidate) {
        replace_candidate(weakest_candidate(), new_candidate);
    }
    void eliminate_candidate(Candidate to_eliminate, bool restore_tree_structure = true) {
        size_t external_index = to_eliminate.winner_tree_external_index;
        if (tree[external_index].candidate != to_eliminate)
            throw std::logic_error("eliminate_candidate: candidate does not belong to this tree");
        tree[external_index].candidate.eliminate();
        //if this tree does not correspond to the active query, theres no reason to maintain tree structure
        //(no more potential incomming candidates, we're just waiting to output the non-eliminated ones)
        if (restore_tree_structure)
            replay_up(external_index);
    }
    void reset() {
        //this function allows us to reuse this class instance for a new query
        //reset the left K external nodes
        iterator iter_end = end();
        for (iterator iter = begin(); iter != iter_end; ++iter)
            iter->reset();
        //bubble the reset nodes up the tree
        prepare_subtree(1);
    }
};

//the following class allows us to detect and mitigate trivial matches in constant time.
//there should only be a single instance, since trivial matches are only relevant within the past QUERY_LEN queries.
class QueryLookup {
private:
    struct CandidateBucket {
        //buckets only have space for a single candidate since we only hold buckets for
        //the most recent QUERY_LEN queries
        Candidate candidate;
        bool contains_candidate;
        bool detect_trivial_match(size_t other_candidate_position) const {
            return contains_candidate && std::abs(candidate.loc - other_candidate_position) < QUERY_LEN;
        }
        void pop_candidate() {
            contains_candidate = false;
        }
        void push_candidate(Candidate c) {
            contains_candidate = true;
            candidate = c;
        }
    };
    CandidateBucket* buckets;
    size_t num_buckets;
    size_t get_bucket_index(size_t candidate_position) const {
        return candidate_position / QUERY_LEN;
    }
    bool detect_trivial_match_bucket(const Candidate in_new, CandidateBucket& out_old) const {
        //todo: verify this logic - would there ever be two candidates in the vacinity?
        size_t new_candidate_position = in_new.loc;
        size_t bucket_index = get_bucket_index(new_candidate_position);
        if (bucket_index > 0 && buckets[bucket_index - 1].detect_trivial_match(new_candidate_position)) {
            out_old = buckets[bucket_index - 1];
            return true;
        }
        if (buckets[bucket_index].detect_trivial_match(new_candidate_position)) {
            out_old = buckets[bucket_index];
            return true;
        }
        if (bucket_index + 1 < num_buckets && buckets[bucket_index + 1].detect_trivial_match(new_candidate_position)) {
            out_old = buckets[bucket_index + 1];
            return true;
        }
        return false;
    }
public:
    QueryLookup() {
        num_buckets = ceil((float)TIME_SERIES_LEN / QUERY_LEN);
        buckets = new CandidateBucket[num_buckets];
    }
    ~QueryLookup() {
        delete[] buckets;
    }
    bool detect_trivial_match(const Candidate in_new, Candidate& out_old) const {
        CandidateBucket bucket;
        if (detect_trivial_match_bucket(in_new, bucket)) {
            out_old = bucket.candidate;
            return true;
        }
        return false;
    }
    //remove deletes all candidates in the vacinity of old_candidate
    void pop_candidate(const Candidate old_candidate) {
        size_t bucket_index = get_bucket_index(old_candidate.loc);
        buckets[bucket_index].pop_candidate();
    }
    void push_candidate(const Candidate new_candidate) {
        //note: this indiscriminately replaces a previous candidate with the new candidate, regardless
        //of which one is a better match (that responsibility is left to the client)
        CandidateBucket old_candidate_bucket;
        if (detect_trivial_match_bucket(new_candidate, old_candidate_bucket))
            old_candidate_bucket.pop_candidate();
        size_t bucket_index = get_bucket_index(new_candidate);
        buckets[bucket_index].push_candidate(new_candidate);
    }
};

class MotifFinder {
private:
    const double* time_series;
    Cache cache;
    //top_k_query_candidates is a circular buffer that holds the best candidates for the last QUERY_LEN queries.
    //it is used to detect and ignore trivial matches.
    TopCandidates* top_k_query_candidates;
    QueryLookup query_lookup;
    std::ofstream results_output;
    const char field_separator;
    size_t curr_query_loc;
    TopCandidates& get_top_k(size_t query_loc) {
        size_t circ_buff_index = query_loc % QUERY_LEN;
        TopCandidates& top_k = query_candidates[circ_buff_index];
        return top_k;
    }
    void print_and_disregard_top_k(size_t query_loc) {
        TopCandidates& top_k = get_top_k(query_loc);
        TopCandidates::iterator iter_end = top_k.end();
        for (auto iter = top_k.begin(); iter != iter_end; ++iter) {
            Candidate c = *iter;
            results_output << c.query_loc << field_separator << c.loc << field_separator << c.dist << std::endl;
            //no more trivial matches to this candidate
            query_lookup.pop_candidate(c);
        }
        top_k.reset();
    }
public:
    MotifFinder(const double* ts, std::string output_filepath, const char results_field_separator = '\t'):
        time_series(ts),
        cache(time_series),
        results_output(output_filepath, std::ofstream::out | std::ofstream::trunc),
        field_separator(results_field_separator)
    {
        query_candidates = new TopCandidates[QUERY_LEN];
    }
    ~MotifFinder() {
        delete[] query_candidates;
    }
    Candidate::DIST_TYPE weakest_distance() {
        TopCandidates top_k = get_top_k(curr_query_loc);
        Candidate weakest_candidate = top_k.weakest_candidate;
        return weakest_candidate.dist;
    }
    void register_candidate(size_t new_candidate_index, Candidate::DIST_TYPE distance) {
        //note: this function does not check to see that the new candidate is actually better and worth remembering
        Candidate new_candidate;
        new_candidate.dist = distance;
        new_candidate.query_loc = curr_query_loc;
        new_candidate.loc = new_candidate_index;
        Candidate old_candidate;
        TopCandidates& curr_top_k = get_top_k(curr_query_loc);
        if (query_lookup.detect_trivial_match(new_candidate, old_candidate)) {
            //found a trivial match
            if (new_candidate.dist < old_candidate.dist) {
                //old candidate is trivial match 
                query_lookup.push_candidate(new_candidate);
                if (old_candidate.query_loc == new_candidate.query_loc) {
                    //trivial match to candidate of current query
                    curr_top_k.replace_candidate(old_candidate, new_candidate);
                } else {
                    //trivial match to candidate of previous query
                    TopCandidates& prev_top_k = get_top_k(old_candidate.query_loc);
                    prev_top_k.eliminate_candidate(old_candidate, false);
                }
            }
        } else {
            //no trivial match, just replace weakest match on current query
            curr_top_k.replace_weakest_candidate(new_candidate);
        }
    }
};

#endif
