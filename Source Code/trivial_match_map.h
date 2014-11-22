#ifndef _TRIVIAL_MATCH_MAP_H_
#define _TRIVIAL_MATCH_MAP_H_

#include "candidate.h"
#include "common.h"

//the following class uses a hash table (hashing function is position/QUERY_LEN) and allows us to detect and mitigate
//trivial matches in constant time. in addition, a linked list is maintained chaining together active hash table items
//so that they can be quickly traversed and deleted, allowing us to reset/reuse the class instance quickly
class TrivialMatchMap {
private:
    struct CandidateBucket {
        Candidate candidate;
        bool contains_candidate;

        //maintain a bucket list for fast reset()
        size_t prev_bucket_index;
        bool points_to_prev_bucket;
        size_t next_bucket_index;
        bool points_to_next_bucket;
        bool detect_trivial_match(Candidate const& other_candidate) const {
            return contains_candidate && candidate.is_candidate_close_to(other_candidate);
        }
        void pop_candidate() {
            contains_candidate = false;
        }
        void push_candidate(Candidate c) {
            contains_candidate = true;
            candidate = c;
        }
        CandidateBucket(): contains_candidate(false), points_to_prev_bucket(false), points_to_next_bucket(false) {}
    };
    CandidateBucket bucket_list_head; //this acts as a dummy node in the bucket list
    const size_t num_buckets;
    std::vector<CandidateBucket> buckets;
    size_t get_bucket_index(size_t candidate_position) const {
        return candidate_position / QUERY_LEN;
    }
    CandidateBucket* get_trivial_match_bucket(Candidate const& new_candidate) {
        //todo: verify this logic - would there ever be two candidates in the vacinity?
        size_t new_candidate_position = new_candidate.loc;
        size_t bucket_index = get_bucket_index(new_candidate_position);
        if (bucket_index > 0 && buckets[bucket_index - 1].detect_trivial_match(new_candidate)) {
            return &buckets[bucket_index - 1];
        }
        if (buckets[bucket_index].detect_trivial_match(new_candidate)) {
            return &buckets[bucket_index];
        }
        if (bucket_index + 1 < num_buckets && buckets[bucket_index + 1].detect_trivial_match(new_candidate)) {
            return &buckets[bucket_index + 1];
        }
        return nullptr;
    }
    void empty_bucket(CandidateBucket* bucket) {
        bucket->pop_candidate();
        if (bucket->points_to_prev_bucket) {
            CandidateBucket& prev_bucket = buckets[bucket->prev_bucket_index];
            prev_bucket.points_to_next_bucket = bucket->points_to_next_bucket;
            prev_bucket.next_bucket_index = bucket->next_bucket_index;
        }
        if (bucket->points_to_next_bucket) {
            CandidateBucket& next_bucket = buckets[bucket->next_bucket_index];
            next_bucket.points_to_prev_bucket = bucket->points_to_prev_bucket;
            next_bucket.prev_bucket_index = bucket->prev_bucket_index;
        }
    }
public:
    TrivialMatchMap():
        num_buckets(ceil(static_cast<float>(TIME_SERIES_LEN) / QUERY_LEN)),
        buckets(num_buckets)
    {}
    bool detect_trivial_match(Candidate const& new_candidate, Candidate* old_candidate) {
        CandidateBucket* bucket = get_trivial_match_bucket(new_candidate);
        if (bucket) {
            *old_candidate = bucket->candidate;
            return true;
        }
        return false;
    }
    void pop_candidate(Candidate const& old_candidate) {
        size_t bucket_index = get_bucket_index(old_candidate.loc);
        empty_bucket(&buckets[bucket_index]);
    }
    void push_candidate(Candidate const& new_candidate) {
        //note: this indiscriminately replaces a previous candidate with the new candidate, regardless
        //of which one is a better match (that responsibility is left to the client)
        CandidateBucket* old_candidate_bucket = get_trivial_match_bucket(new_candidate);
        if (old_candidate_bucket)
            old_candidate_bucket->pop_candidate();
        size_t bucket_index = get_bucket_index(new_candidate.loc);
        CandidateBucket& bucket = buckets[bucket_index];
        bucket.push_candidate(new_candidate);
        if (bucket_list_head.points_to_next_bucket) {
            bucket.next_bucket_index = bucket_list_head.next_bucket_index;
            bucket.points_to_next_bucket = true;
            CandidateBucket& other_bucket = buckets[bucket_list_head.next_bucket_index];
            other_bucket.prev_bucket_index = bucket_index;
            other_bucket.points_to_prev_bucket = true;
        }
        bucket_list_head.next_bucket_index = bucket_index;
        bucket_list_head.points_to_next_bucket = true;
    }
    void reset() {
        //delete all buckets by traversing the list containing the active ones.
        //bucket_list_head is not contained within the bucket vector, so it is *not* pointed to
        //by the first node in the bucket list. therefore, we use a pointer to iterate over all
        //the buckets and delete them.
        CandidateBucket* bucket_ptr = &bucket_list_head;
        while (bucket_ptr->points_to_next_bucket) {
            bucket_ptr = &buckets[bucket_ptr->next_bucket_index];
            empty_bucket(bucket_ptr);
        }
        bucket_list_head.points_to_next_bucket = false;
    }
};

#endif
