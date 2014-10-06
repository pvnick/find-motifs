//config
#define USE_MPI

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"

#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iterator>
#include <cmath>
#include <boost/mpi.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/asio.hpp>
#include <cstdio>
#include <stdlib.h>
#include <ctime>
#include <csignal>
#include <memory>
#include <vector>


namespace mpi = boost::mpi;
namespace interprocess = boost::interprocess;
#ifdef USE_MPI
    static bool use_mpi = true;
#else
    static bool use_mpi = false;
#endif

unsigned int query_start_pos() {
    if (use_mpi) {
        mpi::communicator world;
        return world.rank();
    } else {
        return 0;
    }
}

unsigned int search_increment() {
    if (use_mpi) {
        mpi::communicator world;
        return world.size();
    } else {
        return 1;
    }
}

std::string get_output_filename() {
    if (use_mpi) {
        mpi::communicator world;
        std::ostringstream filename;
        filename << "results" << world.rank() << ".out";
        return filename.str();
    } else {
        return "results.out";
    }
}

std::ostream& msg(std::string str) {
    if (use_mpi) {
        mpi::communicator world;
        std::cerr << "Proc " << world.rank() << ": " << str;
    } else {
        std::cerr << str;
    }
    return std::cerr;
}

std::ostream& msgl(std::string str) {
    return msg(str) << std::endl;
}

#include "lemire_envelope.h"
#include "cache.h"
#include "find_motifs.h"

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
