#include "common.h"

#include "find_motifs.h"
#include <cstdio>
#include <csignal>

#ifdef USE_PROFILER
    #include "profiler.h"
#endif

void sighandler(int sig)
{
    //dont leave shared memory regions lying around (I don't know
    //what would happen, but it doesn't sound like a good idea)
    cache_type::destroy();
    exit(1);
}

/// Main Function
int main(int argc, char *argv[])
{
#ifdef USE_MPI
    mpi::environment env(argc, argv);
#endif
    unsigned int K = 100;

    size_t start_pos = query_start_pos();
    size_t increment = search_increment();
#ifdef USE_PROFILER
    ProfilerStart("/tmp/profile");
#endif

    std::cerr << "Reading time series data file" << std::endl;
    double* time_series;
    std::ifstream in(SERIES_FILEPATH);
    time_series = new double[TIME_SERIES_LEN];
    double point;
    for (int i = 0; in >> point && i != TIME_SERIES_LEN; ++i)
        time_series[i] = point;
    in.close();
    cache_type::init(time_series);

    signal(SIGABRT, &sighandler);
	signal(SIGTERM, &sighandler);
	signal(SIGINT, &sighandler);

    MotifFinder engine;

    std::ofstream out(std::string("/scratch/lfs/pvnick/motif_results/") + get_output_filename(),
                      std::ofstream::out | std::ofstream::trunc);
    for (size_t i = start_pos
         ; i < TIME_SERIES_LEN - QUERY_LEN /*don't query at the end of the time series*/ //todo: should this be 2 x QUERY_LEN?
         ; i += increment)
    {
        msg("Querying from ") << i << std::endl;
        MotifFinder::TopKMatches results = engine.single_pass(K, i);
        results >> out;
        msgl("Query complete");
    }
    out.close();
#ifdef USE_PROFILER
    ProfilerStop();
#endif
    return 0;
}
