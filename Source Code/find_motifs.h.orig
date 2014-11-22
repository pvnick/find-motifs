#ifndef _FIND_MOTIFS_H_
#define _FIND_MOTIFS_H_

#include "common.h"
#include "cache.h"
#include "ucr_dtw.h"
#include <limits>

/*note to self:
call signatures should generally be one of the following:
bar(foo f); // want to obtain a copy of f
bar(const foo& f); // want to read f
bar(foo* f); // want to modify f
*/

class TopCandidates;

struct Candidate {
private:
    static constexpr double max_dist = std::numeric_limits<dist_type>::max();
public:
    dist_type dist;
    size_t loc;
    size_t query_loc;
    size_t winner_tree_external_index; //holding this allows faster updating of the top-k tree
    //if the candidate held in this node is, at a subsequent query position, determined to be a trivial match,
    //is_disabled will get set to true, and this candidate will be skipped when candidates are outputted
    bool is_disabled;
    Candidate() {
        reset();
    };
    void disable() {
        reset();
        is_disabled = true;
    }
    bool is_valid() const {
        //since candidates must be at least QUERY_LEN points past the start of the time series, their location can't be zero
        return loc != 0;
    }
    void reset() {
        is_disabled = false;
        dist = max_dist;
        loc = 0;
        query_loc = 0;
    }
    bool operator==(const Candidate& rhs) {
        return query_loc == rhs.query_loc && loc == rhs.loc;
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
        bool operator==(const Node& rhs) {
            return is_placeholder == rhs.is_placeholder || candidate == rhs.candidate;
        }
    };
    //note: packed array is 1-based to make the math easier
    std::vector<Node> tree;
    size_t num_nodes;
    size_t K;
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
        const Node& left_child = tree[2 * index];
        const Node& right_child = tree[2 * index + 1];
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
    size_t count_internal_nodes() const {
        return (num_nodes - 1) / 2;
    }
public:
    class TopCandidates_Iter: public std::iterator<std::forward_iterator_tag, Candidate> {
    //we don't care about order when iterating through the top-k matches, so just traverse the winner tree external nodes
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
        size_t first_external_placeholder = first_external_index + K;
        iterator iter(&tree, first_external_placeholder);
        return iter;
    }
    TopCandidates(size_t max_entries): K(max_entries) {
        //we need 2^n external nodes for easy math, where n is an integer
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
    }
    Candidate weakest_candidate() const {
        return tree[1].candidate;
    }
    size_t weakest_candidate_external_index() const {
        const Candidate& m = weakest_candidate();
        return m.winner_tree_external_index;
    }
    void replace_candidate(Candidate const& old_candidate, Candidate const& new_candidate) {
        size_t external_index = old_candidate.winner_tree_external_index;
        tree[external_index].candidate = new_candidate;
        tree[external_index].candidate.winner_tree_external_index = external_index;
        replay_up(external_index);
    }
    void replace_weakest_candidate(Candidate const& new_candidate) {
        replace_candidate(weakest_candidate(), new_candidate);
    }
    void disable_candidate_external_node(Candidate& to_disable) {
        //if this tree does not correspond to the active query - only time this function will be called - theres
        //no reason to maintain tree structure (no more potential incomming candidates, we're just waiting to
        //output the non-disabled ones) so this function should not be called if this class instance
        //belongs to the current query; replace_candidate should be used instead
        size_t external_index = to_disable.winner_tree_external_index;
        assert(tree[external_index].candidate == to_disable);
        tree[external_index].candidate.disable();
    }
    void reset() {
        //this function allows us to reuse this class instance for a new query.
        //reset the left K external nodes
        iterator iter_end = end();
        for (iterator iter = begin(); iter != iter_end; ++iter)
            iter->reset();
        //bubble the reset nodes up the tree
        prepare_subtree(1);
    }
};

//the following class uses a hash table (hashing function is position/QUERY_LEN) and allows us to detect and mitigate
//trivial matches in constant time. there should only be a single instance, since trivial matches are only relevant
//within the past QUERY_LEN queries.
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
    std::vector<CandidateBucket> buckets;
    size_t num_buckets;
    size_t get_bucket_index(size_t candidate_position) const {
        return candidate_position / QUERY_LEN;
    }
    CandidateBucket* get_trivial_match_bucket(Candidate const& new_candidate) {
        //todo: verify this logic - would there ever be two candidates in the vacinity?
        size_t new_candidate_position = new_candidate.loc;
        size_t bucket_index = get_bucket_index(new_candidate_position);
        if (bucket_index > 0 && buckets[bucket_index - 1].detect_trivial_match(new_candidate_position)) {
            return &buckets[bucket_index - 1];
        }
        if (buckets[bucket_index].detect_trivial_match(new_candidate_position)) {
            return &buckets[bucket_index];
        }
        if (bucket_index + 1 < num_buckets && buckets[bucket_index + 1].detect_trivial_match(new_candidate_position)) {
            return &buckets[bucket_index + 1];
        }
        return nullptr;
    }
public:
    QueryLookup() {
        num_buckets = ceil(static_cast<float>(TIME_SERIES_LEN) / QUERY_LEN);
        buckets.resize(num_buckets);
    }
    bool detect_trivial_match(const Candidate& new_candidate, Candidate* old_candidate) {
        CandidateBucket* bucket = get_trivial_match_bucket(new_candidate);
        if (bucket) {
            *old_candidate = bucket->candidate;
            return true;
        }
        return false;
    }
    void pop_candidate(Candidate const& old_candidate) {
        size_t bucket_index = get_bucket_index(old_candidate.loc);
        buckets[bucket_index].pop_candidate();
    }
    void push_candidate(Candidate const& new_candidate) {
        //note: this indiscriminately replaces a previous candidate with the new candidate, regardless
        //of which one is a better match (that responsibility is left to the client)
        CandidateBucket* old_candidate_bucket = get_trivial_match_bucket(new_candidate);
        if (old_candidate_bucket)
            old_candidate_bucket->pop_candidate();
        size_t bucket_index = get_bucket_index(new_candidate.loc);
        buckets[bucket_index].push_candidate(new_candidate);
    }
};

class MotifFinder {
private:
    std::vector<double> time_series;
    QueryLookup query_lookup;
    std::ofstream results_output;
    const char field_separator;
    //top_k_query_candidates is a circular buffer that holds the best candidates for the last QUERY_LEN queries.
    //it is used to detect and ignore trivial matches.
    std::vector<TopCandidates> top_k_query_candidates;
    Cache cache;
    size_t curr_query_loc;
    TopCandidates& get_top_k(size_t query_loc) {
        size_t circ_buff_index = query_loc % QUERY_LEN;
        TopCandidates& top_k = top_k_query_candidates[circ_buff_index];
        return top_k;
    }
    void output_and_disregard_top_k(size_t query_loc) {
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
    MotifFinder(std::string const& series_filepath, size_t series_len, unsigned int K, std::string const& output_filepath, const char results_field_separator):
        time_series(series_len),
        results_output(output_filepath, std::ofstream::out | std::ofstream::trunc),
        field_separator(results_field_separator),
        top_k_query_candidates(QUERY_LEN, TopCandidates(K)),
        cache(time_series),
        curr_query_loc(0)
    {
        std::cerr << "Reading time series data file" << std::endl;
        std::ifstream in(series_filepath);
        double point;
        while (in >> point)
            time_series.push_back(point);
        in.close();
    }
    dist_type weakest_distance() {
        TopCandidates& top_k = get_top_k(curr_query_loc);
        const Candidate& weakest_candidate = top_k.weakest_candidate();
        return weakest_candidate.dist;
    }
    void try_register_candidate(size_t new_candidate_index, dist_type distance) {
        //note: this function does not check to see that the new candidate is actually better and worth remembering.
        //that responsibility is left to the client, who should be checking that with weakest_distance
        TopCandidates& curr_top_k = get_top_k(curr_query_loc);
        Candidate new_candidate;
        new_candidate.dist = distance;
        new_candidate.query_loc = curr_query_loc;
        new_candidate.loc = new_candidate_index;
        Candidate old_candidate;
        if (query_lookup.detect_trivial_match(new_candidate, &old_candidate)) {
            //found a trivial match
            if (new_candidate.dist < old_candidate.dist) {
                //old candidate is trivial match (new is better), so disregard old and register the new one
                query_lookup.push_candidate(new_candidate);
                if (old_candidate.query_loc == new_candidate.query_loc) {
                    //trivial match to candidate of current query
                    curr_top_k.replace_candidate(old_candidate, new_candidate);
                } else {
                    //trivial match to candidate of previous query
                    TopCandidates& prev_top_k = get_top_k(old_candidate.query_loc);
                    prev_top_k.disable_candidate_external_node(old_candidate);
                }
            }
        } else {
            //no trivial match
            query_lookup.push_candidate(new_candidate);
            curr_top_k.replace_weakest_candidate(new_candidate);
        }
    }
    void run(size_t start_pos, size_t end_pos, size_t candidate_increment) {
        using namespace std::placeholders;
        UCR_DTW subsequence_search(&cache);
        auto weakest_dist_callback = std::bind(&MotifFinder::weakest_distance, this);
        auto register_candidate_callback = std::bind(&MotifFinder::try_register_candidate, this, _1, _2);
        for (size_t i = start_pos; i != end_pos; ++i) {
            curr_query_loc = i;
            subsequence_search.single_pass(i, candidate_increment, weakest_dist_callback, register_candidate_callback);
            if (i - start_pos >= QUERY_LEN) {
                //no more trivial to oldest query stored in the circular buffer
                output_and_disregard_top_k(i - QUERY_LEN);
            }
        }
/*
 32 void do_call(std::function<void(std::string, std::string)> f) {
 33     f("hello globe", "goodbye");
 34 }
 35
 36 void callfunc() {
 37     using namespace std::placeholders;  // for _1, _2, _3...
 38     foo f;
 39     auto f3 = std::bind(&foo::testfunc, &f, _1, _2);
 40     do_call(f3);
 41 }
 42
*/
    }
};

#endif
