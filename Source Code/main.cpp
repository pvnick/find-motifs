#include "find_motifs.h"
#include <iostream>
#include <string>

size_t query_start_position(size_t series_length, unsigned int proc_rank, unsigned int num_procs) {
    size_t search_space_length = floor((float)proc_rank / num_procs * series_length * (series_length + 1.0) / 2.0);
    size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
    return position;
}

/// Main Function
int main(  int argc , char *argv[] )
{
    for (int i = 0; i != 8; ++i) {
        size_t start = query_start_position(TIME_SERIES_LEN, i, 8);
        size_t fin = query_start_position(TIME_SERIES_LEN, i + 1, 8);
        std::cout << (TIME_SERIES_LEN - fin) * (TIME_SERIES_LEN - fin + 1) / 2 - (TIME_SERIES_LEN - start) * (TIME_SERIES_LEN - start + 1) / 2 << std::endl;
    }
    exit(0);
    unsigned int K = 100;
    MotifFinder engine;
    for (int i = 0; i != 100; ++i) {
        MotifFinder::TopKMatches result = engine.single_pass(K, i);
        result.print();
    }
    return 0;
}
