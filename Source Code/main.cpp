//#define USE_MPI
//#define USE_PROFILER

#include "find_motifs.h"
#include <iostream>
#include <string>

#ifdef USE_PROFILER
    #include "profiler.h"
#endif
#ifdef USE_MPI
    #include <boost/mpi.hpp>
    #include <boost/serialization/string.hpp>

    namespace mpi = boost::mpi;

    size_t query_start_position_by_params(size_t series_length, unsigned int proc_rank, unsigned int num_procs) {
        proc_rank = num_procs - proc_rank; //make assignment easier to conceptualize
        size_t search_space_length = floor((float)proc_rank / num_procs * series_length * (series_length + 1.0) / 2.0);
        size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
        return position;
    }

    size_t my_query_start_pos() {
        mpi::environment env;
        mpi::communicator world;
        unsigned int num_procs = world.size();
        unsigned int proc_rank = world.rank();
        return query_start_position_by_params(TIME_SERIES_LEN, proc_rank, num_procs);
    }

    size_t my_query_end_pos() {
        mpi::environment env;
        mpi::communicator world;
        unsigned int num_procs = world.size();
        unsigned int proc_rank = world.rank();
        return query_start_position_by_params(TIME_SERIES_LEN, proc_rank + 1, num_procs);
    }
#else
    size_t my_query_start_pos() {
        return 0;
    }

    size_t my_query_end_pos() {
        return TIME_SERIES_LEN;
    }
#endif

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
    */
    unsigned int K = 100;
    size_t start_pos = my_query_start_pos();
    size_t end_pos = my_query_end_pos();
    MotifFinder engine;
    #ifdef USE_PROFILER
        ProfilerStart("/tmp/profile");
        end_pos = start_pos + 3;
    #endif
    for (size_t i = start_pos; i != end_pos
                            && i != TIME_SERIES_LEN - QUERY_LEN /*don't query at the end of the time series*/
                            ; ++i) {
        std::cout << "Querying from " << i << std::endl;
        MotifFinder::TopKMatches result = engine.single_pass(K, i);
        result.print();
    }
    #ifdef USE_PROFILER
        ProfilerStop();
    #endif
    return 0;
}
