#define USE_MPI
//#define USE_PROFILER

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"


#include "find_motifs.h"
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <stdexcept>

#ifdef USE_PROFILER
    #include "profiler.h"
#endif

/// Main Function
int main(int argc, char *argv[])
{
#ifdef USE_MPI
    mpi::environment env(argc, argv);
    construct_hostname_proc_layouts();
#endif
    unsigned int K = 100;
    size_t start_pos = query_start_pos();
    size_t increment = search_increment();
    std::shared_ptr<MotifFinder> engine = get_engine();
#ifdef USE_PROFILER
    ProfilerStart("/tmp/profile");
#endif
    std::ofstream out(std::string("/scratch/lfs/pvnick/motif_results/") + get_output_filename(),
                      std::ofstream::out | std::ofstream::trunc);
    for (size_t i = start_pos;
         i < TIME_SERIES_LEN - QUERY_LEN /*don't query at the end of the time series*/; //todo: should this be 2 x QUERY_LEN?
         i += increment)
    {
        msg("Querying from ") << i << std::endl;
        MotifFinder::TopKMatches results = engine->single_pass(K, i);
        results >> out;
        msgl("Query complete");
    }
    out.close();
#ifdef USE_PROFILER
    ProfilerStop();
#endif
    return 0;
}
