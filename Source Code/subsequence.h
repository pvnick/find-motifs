//todo: I think this is possible to do without shared memory!
//
#ifndef _SUBSEQUENCE_H_
#define _SUBSEQUENCE_H_

#include "common.h"
#include "lemire_envelope.h"
#include <limits>

class Subsequence {
public:
    size_t time_series_pos;
    std::vector<double> series_normalized;
    double mean;
    double stddev;
    double range;
    size_t fragment_length;
    LemireEnvelope lemire_envelope;
    bool initialized;

    void init(const std::vector<double>& time_series, size_t position) {
        time_series_pos = position;
        mean = 0;
        stddev = 0;
        series_normalized.resize(QUERY_LEN);
        double ex = 0, ex2 = 0;
        double min = std::numeric_limits<double>::max(),
               max = std::numeric_limits<double>::min();
        for (fragment_length = 0; fragment_length + time_series_pos != TIME_SERIES_LEN && fragment_length != QUERY_LEN; ++fragment_length) {
            double d = time_series[fragment_length + time_series_pos];
            ex += d;
            ex2 += d*d;
        }
        mean = ex/fragment_length;
        stddev = ex2/fragment_length;
        stddev = sqrt(stddev-mean*mean);
        for(size_t i = 0; i != fragment_length; i++) {
             double d = time_series[i + time_series_pos];
             series_normalized[i] = (d - mean)/stddev;
             if (d > max) max = d;
             if (d < min) min = d;
        }
        lemire_envelope = LemireEnvelope(time_series, time_series_pos, WARPING_r);
        range = max - min;
        initialized = true;
    }

    Subsequence(): initialized(false) {}
};

class SubsequenceLookup {
private:
    std::vector<Subsequence> entries;
    const std::vector<double>& time_series;
public:
    SubsequenceLookup() = delete;
    SubsequenceLookup(const std::vector<double>& ts): time_series(ts) {
        entries.resize(ts.size());
    }
    Subsequence const& operator[](size_t position) {
        //memoization lets us compute only relevant subsequence entries and then reuse them
        if ( ! entries[position].initialized)
            entries[position].init(time_series, position);
        return entries[position];
    }
    Subsequence operator=(Subsequence const& rhs) = delete;
};

#endif // _SUBSEQUENCE_H_
