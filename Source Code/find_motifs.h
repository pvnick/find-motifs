#ifndef _FIND_MOTIFS_H_
#define _FIND_MOTIFS_H_

#include "common.h"
#include "cache.h"
#include "ucr_dtw.h"
#include <limits>

struct Candidate {
public:
    typedef double DIST_TYPE;
private:
    static constexpr double max_dist = std::numeric_limits<DIST_TYPE>::max();
public:
    DIST_TYPE dist;
    size_t loc;
    size_t winner_tree_external_index;
    Candidate(): dist(max_dist) {};
};

//the following class uses a winner tree to hold the top-K items
class TopCandidates {
private:
    struct Node {
        Candidate candidate;
        bool is_placeholder; //indicates whether this holds a real candidate or is just an array slot placeholder
    };
    //note: packed array is 1-based to make the math easier (ie slot 0 is unused)
    Node* tree;
    size_t num_nodes;
    const bool is_internal(size_t index) const {
        return index <= num_nodes / 2;
    }
    size_t prepare_tree(size_t index) {
        //since all slots represent candidates which are initially the same distance, we just need to
        //bubble up stored external indices (shaves lg(N) comparisons off replacing the weakest candidate)
        //and whether the node is a placeholder. this function is recursive and returns an index to an external node slot.
        //we use the left leg for this, so this function returns the index to the external node arrived
        //at if we were to follow the left leg from the current position to the bottom of the tree.
        //consequentially, adding items preserves tree completeness, and we don't need manually ensure that
        //we're storing no more than K items
        size_t external_index;
        size_t left_child_index = index * 2;
        size_t right_child_index = left_child_index + 1;
        if (is_internal(left_child_index)) {
            //if left child is an internal node, then right child is guaranteed to be internal
            external_index = prepare_tree(left_child_index);
            prepare_tree(right_child_index);
        } else {
            //this is an external node, so return the current index, which may or may not be relevant to the caller
            external_index = index;
        }
        tree[index] = tree[external_index]; //this will bubble placeholder designations up to where the left sibling is a non-placeholder
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
public:
    TopCandidates(unsigned int K) {
        //we need 2^n external nodes for easy math, where n is an integer
        size_t num_external_nodes = pow(2, ceil(lg(K)));
        size_t num_internal_nodes = num_external_nodes - 1;
        num_nodes = num_external_nodes + num_internal_nodes;
        tree = new Node[num_nodes];
        //since we can only store K candidates in the tree, disable external nodes more than K slots from the left.
        //doing so automagically guarantees we store at most K items in the tree
        for (int i = num_nodes; i > num_internal_nodes + K; --i)
            tree[i].is_placeholder = true;
        prepare_tree(1);
    }
    ~TopCandidates() {
        delete[] tree;
    }
    const Candidate weakest_candidate() const {
        return tree[1].candidate;
    }
    const size_t weakest_candidate_external_index() const {
        Candidate m = weakest_candidate();
        return m.winner_tree_external_index;
    }
    void replay_up(size_t external_index) {
        for(size_t i = external_index / 2; i != 0; i /= 2) {
            Node winner = get_winning_child(i);
            if (tree[i] == winner)
                break; //contest outcome didn't change
            tree[i] = winner;
        }
    }
    void replace_weakest_candidate(Candidate replacement) {
        size_t external_index = winner_external_node_index();
        replace_candidate(external_index, replacement);
    }
    void replace_candidate(size_t external_index, Candidate replacement) {
        replacement.winner_tree_external_index = external_index;
        tree[external_index].candidate = replacement;
        replay_up(external_index);
    }
};

//the following class allows us to detect and mitigate trivial matches in constant time
class QueryLookup {
private:
    struct CandidateBucket {
        //buckets only have space for a single candidate since we only hold buckets for
        //the most recent QUERY_LEN queries
        Candidate* candidate;
        CandidateBucket: candidate(nullptr) {}
    };
    CandidateBucket* buckets;
    size_t num_buckets;
public:
    QueryLookup() {
        num_buckets = ceil((float)TIME_SERIES_LEN / QUERY_LEN);
        buckets = new CandidateBucket[num_buckets];
    }
    ~QueryLookup() {
        delete[] buckets;
    }
    Candidate* lookup_previous_candidate(size_t candidate_position) {
        size_t bucket_index = candidate_position / QUERY_LEN;
        CandidateBucket bucket = buckets[bucket_index];
        //right here, should be checking more than 1 bucket slot
        return bucket.candidate;
    }
    void replace_candidate(size_t candidate_position) {

    }
};

/*
    const size_t query_pos;
    const size_t max_size;
    size_t curr_size;
    const size_t query_length;
    Match* matches_head;
    void delete_weakest_match() {
        Match* m = matches_head->next;
        matches_head->next = m->next;
        delete m;
        --curr_size;
    }
    void shrink_if_necessary() {
        if (size() > max_size) {
            delete_weakest_match();
        }
    }
public:
    size_t size() {
        return curr_size;
    }
    double weakest_dist() {
        if (size() < max_size) return INF;
        return matches_head->next->dist;
    }
    TopKMatches(const size_t k, const size_t len, const size_t pos):
        query_pos(pos),
        max_size(k),
        curr_size(0),
        query_length(len)
    {
        matches_head = new Match();
        matches_head->is_dummy = true;
        matches_head->next = nullptr;
    };
    ~TopKMatches() {
        clear();
        delete matches_head;
    }
    TopKMatches(const TopKMatches& src): TopKMatches(src.max_size, src.query_length, src.query_pos) {
        for (Match* m = src.matches_head->next; m != nullptr; m = m->next)
            insert_match(m->loc, m->dist);
    }
    void clear() {
        while (size()) {
            delete_weakest_match();
        }
    }
    Match* get_self_match_preceeding(long long loc) {
        Match* m = matches_head;
        for (m = matches_head; m->next != nullptr && (size_t)abs(m->next->loc - loc) >= query_length; m = m->next);
        return m;
    }
    void insert_match(long long loc, double dist) {

        Match* self_match_preceeding = get_self_match_preceeding(loc);
        if (self_match_preceeding->next != nullptr) {
            if (dist < self_match_preceeding->next->dist) {
                //delete the self-match, then insert in to the correct slot
                Match* self_match = self_match_preceeding->next;
                self_match_preceeding->next = self_match->next;
                delete self_match;
                --curr_size;
            } else {
                //found self-match which is weaker than what we already got, so ignore it
                return;
            }
        }

        Match* preceeding = matches_head;
        for (; preceeding->next != nullptr && preceeding->next->dist > dist; preceeding = preceeding->next);
        Match* new_match = new Match();
        new_match->loc = loc;
        new_match->dist = dist;
        new_match->next = preceeding->next;
        preceeding->next = new_match;
        ++curr_size;
        shrink_if_necessary();
    }
    std::ostream& operator>>(std::ostream& out) {
        Match* m;
        for (m = matches_head->next; m != nullptr; m = m->next) {
            out << query_pos << "," << m->loc << "," << m->dist << std::endl;
        }
        return out;
    }*/
};


class MotifFinder {
private:
    class TopMatches;
    class Match;
    const double* time_series;
    Cache cache;
    //query_matches is a circular buffer that holds matches for the last QUERY_LEN queries.
    //it is used to detect and ignore trivial matches.
    TopMatches* query_matches;

    size_t candidate_bucket_index(size_t candidate_pos) {
        return candidate_pos / QUERY_LEN;
    }
public:

    MotifFinder(const double* ts):
        time_series(ts),
        cache(time_series)
    {
        query_matches = new TopKMatches[QUERY_LEN];
        size_t num_candidate_buckets = TIME_SERIES_LEN / QUERY_LEN;
        candidate_query_lookup = new Match[num_candidate_buckets];
    }

    ~MotifFinder() {
        delete[] query_matches;
        delete[] candidate_query_lookup;
    }


};

#endif
