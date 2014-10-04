#include "find_motifs.h"
#include <iostream>
#include <string>

size_t query_start_position(size_t series_length, unsigned int proc_rank, unsigned int num_procs) {
    proc_rank = num_procs - proc_rank; //make assignment easier to conceptualize
    size_t search_space_length = floor((float)proc_rank / num_procs * series_length * (series_length + 1.0) / 2.0);
    size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
    return position;
}

/// Main Function
int main(  int argc , char *argv[] )
{
/*
    int num_procs = 100;
    for (int i = 0; i != num_procs; ++i) {
        size_t fin = query_start_position(TIME_SERIES_LEN, i + 1, num_procs);
        size_t start = query_start_position(TIME_SERIES_LEN, i, num_procs);
        std::cout << "Start: " << start << ", End: " << fin << ", Queries: " << fin - start << ", Candidates: ";
        std::cout << (TIME_SERIES_LEN - start) * (TIME_SERIES_LEN - start + 1) / 2 - (TIME_SERIES_LEN - fin) * (TIME_SERIES_LEN - fin + 1) / 2 << std::endl;
    }
    exit(0);*/

    int proc_rank = 0, num_procs = 100;
    unsigned int K = 100;
    size_t query_len = QUERY_LEN;
    //don't query from very end
    size_t search_space = TIME_SERIES_LEN - query_len;
    MotifFinder engine;
    size_t proc_start = query_start_position(search_space, proc_rank, num_procs);
    size_t proc_fin = query_start_position(search_space, proc_rank + 1, num_procs);
    for (int i = proc_start; i <= proc_fin; ++i) {
        std::cout << "Querying from " << i << std::endl;
        MotifFinder::TopKMatches result = engine.single_pass(K, i);
        result.print();
    }

    return 0;
}
