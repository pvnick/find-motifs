#ifndef _CANDIDATE_H_
#define _CANDIDATE_H_

#include "common.h"

struct Candidate {
private:
    static constexpr double max_dist = std::numeric_limits<double>::max();
public:
    double dist;
    size_t loc;
    size_t query_loc;
    size_t winner_tree_external_index; //holding this allows faster updating of the top-k tree
    Candidate() {
        reset();
    };
    bool is_valid() const {
        //since candidates must be at least QUERY_LEN points past the start of the time series, their location can't be zero
        return loc != 0;
    }
    void reset() {
        dist = max_dist;
        loc = 0;
        query_loc = 0;
    }
    bool operator==(const Candidate& rhs) {
        return query_loc == rhs.query_loc && loc == rhs.loc;
    }
    bool is_query_close_to(const Candidate& other) const {
        //todo: verify this is correct and not off by one error
        size_t max_loc = std::max(query_loc, other.query_loc);
        size_t min_loc = std::min(query_loc, other.query_loc);
        return max_loc - min_loc < QUERY_LEN;
    }
    bool is_candidate_close_to(const Candidate& other) const {
        //todo: verify this is correct and not off by one error
        size_t max_loc = std::max(loc, other.loc);
        size_t min_loc = std::min(loc, other.loc);
        return max_loc - min_loc < QUERY_LEN;
    }
};

#endif // _CANDIDATE_H_
