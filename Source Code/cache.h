//todo: I think this is possible to do without shared memory!
//
#ifndef _CACHE_H_
#define _CACHE_H_

#include "common.h"
#include "lemire_envelope.h"
#include <limits>

class CacheEntry {
public:
    unsigned int time_series_pos;
    double series_normalized[QUERY_LEN];
    double mean;
    double stddev;
    double range;
    unsigned int fragment_length;
    LemireEnvelope lemire_envelope;

    CacheEntry(double* time_series, unsigned int position):
        time_series_pos(position),
        mean(0),
        stddev(0)
    {
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
        for(unsigned int i = 0; i != fragment_length; i++) {
             double d = time_series[i + time_series_pos];
             series_normalized[i] = (d - mean)/stddev;
             if (d > max) max = d;
             if (d < min) min = d;
        }
        lemire_envelope = LemireEnvelope(time_series + time_series_pos, WARPING_r);
        range = max - min;
    }
    ~CacheEntry() {}
};

class Cache {
private:
    CacheEntry** entries;
    const double* time_series;
public:
    Cache() = delete;
    Cache(double* ts): time_series(ts) {
        std::cerr << "Cache: Allocating " << TIME_SERIES_LEN * sizeof(CacheEntry) << " bytes of memory for cache entries" << std::endl;
        entries = new CacheEntry*[TIME_SERIES_LEN];
        for (int i = 0; i != TIME_SERIES_LEN; ++i)
            entries[i] = nullptr; //initialize all entries to sentinel value nullptr
    }
    ~Cache() {
        for (int i = 0; i != TIME_SERIES_LEN; ++i)
            if (entries[i] != nullptr)
                delete entries[i];
        delete[] entries;
    }
    CacheEntry& operator[](size_t position) {
        //memoization lets us only compute relevant cache entries and then reuse them
        if (entries[position] == nullptr)
            entries[position] = new Cache(time_series);
        return *(entries[position]);
    }
};

#endif // _CACHE_H_
