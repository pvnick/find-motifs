//config
#define USE_MPI

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"

#include "common.h"
#include "cache.h"
#include "find_motifs.h"

static std::string global_shm_tag;

#ifdef USE_MPI
    bool use_mpi = true;
#else
    bool use_mpi = false;
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

#ifdef USE_PROFILER
    #include "profiler.h"
#endif

static Cache* cache;

//todo: test that this function is even necessary
void sighandler(int sig)
{
    //dont leave shared memory regions lying around (I don't know
    //what would happen, but it doesn't sound like a good idea)
    delete cache;
    exit(1);
}


/// Main Function
int main(int argc, char *argv[])
{
    srand(time(0));
    global_shm_tag = std::to_string((unsigned long long)rand());
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

    signal(SIGABRT, &sighandler);
	signal(SIGTERM, &sighandler);
	signal(SIGINT, &sighandler);

    cache = new cache_type(time_series);
    MotifFinder engine(time_series, *cache);

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
    delete[] time_series;

    out.close();
    delete cache;
#ifdef USE_PROFILER
    ProfilerStop();
#endif
    return 0;
}
