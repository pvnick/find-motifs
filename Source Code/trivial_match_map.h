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
        CandidateBucket* prev_bucket = nullptr;
        CandidateBucket* next_bucket = nullptr;
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
        CandidateBucket(): contains_candidate(false) {}
    };
    CandidateBucket* bucket_list_head = nullptr;
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
        if (bucket->prev_bucket) {
            bucket->prev_bucket->next_bucket = bucket->next_bucket;
        }
        if (bucket->next_bucket) {
            bucket->next_bucket->prev_bucket = bucket->prev_bucket;
        }
        if (bucket_list_head == bucket)
            bucket_list_head = bucket->next_bucket;
        bucket->prev_bucket = bucket->next_bucket = nullptr;
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
        size_t bucket_index = get_bucket_index(new_candidate.loc);
        CandidateBucket* new_bucket = &buckets[bucket_index];
        if (old_candidate_bucket) {
            old_candidate_bucket->pop_candidate();
            //replace the old bucket's list node with the new bucket
            if (old_candidate_bucket->prev_bucket) {
                new_bucket->prev_bucket = old_candidate_bucket->prev_bucket;
                old_candidate_bucket->prev_bucket->next_bucket = new_bucket;
            }
            if (old_candidate_bucket->next_bucket) {
                new_bucket->next_bucket = old_candidate_bucket->next_bucket;
                old_candidate_bucket->next_bucket->prev_bucket = new_bucket;
            }
            if (old_candidate_bucket != new_bucket)
                old_candidate_bucket->prev_bucket = old_candidate_bucket->next_bucket = nullptr;
            if (bucket_list_head == old_candidate_bucket)
                bucket_list_head = new_bucket;
        } else {
            //no trivial match, so no old bucket exists in the linked list.
            //therefore, the new bucket is a *new* node in the linked list
            if (bucket_list_head == nullptr) {
                bucket_list_head = new_bucket;
            } else {
                bucket_list_head->prev_bucket = new_bucket;
                new_bucket->next_bucket = bucket_list_head;
                bucket_list_head = new_bucket;
            }
        }
        new_bucket->push_candidate(new_candidate);
    }
    void reset() {
        while (bucket_list_head) {
            empty_bucket(bucket_list_head);
        }
    }
};

#endif
