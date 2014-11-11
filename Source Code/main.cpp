//config
#define USE_MPI

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define MIN_RANGE 5
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"

#include "common.h"
#include "find_motifs.h"
#include "ucr_dtw.h"

#ifdef USE_MPI
    bool use_mpi = true;
#else
    bool use_mpi = false;
#endif

size_t distributed_query_start_position(size_t series_length, unsigned int proc_rank, unsigned int num_procs) {
    size_t search_space_length = floor((float)(num_procs - proc_rank) / num_procs * series_length * (series_length + 1.0) / 2.0);
    size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
    return position;
}

unsigned int query_start_pos() {
    if (use_mpi) {
        mpi::communicator world;
        return distributed_query_start_position(TIME_SERIES_LEN, world.rank(), world.size());
    } else {
        return 0;
    }
}

unsigned int query_end_pos() {
    if (use_mpi) {
        mpi::communicator world;
        return distributed_query_start_position(TIME_SERIES_LEN, world.rank() + 1, world.size());
    } else {
        return TIME_SERIES_LEN - QUERY_LEN;
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


/// Main Function
int main(int argc, char *argv[])
{
#ifdef USE_MPI
    mpi::environment env(argc, argv);
#endif
    unsigned int K = 100;

    size_t start_pos = query_start_pos();
    size_t end_pos = query_end_pos();
#ifdef USE_PROFILER
    ProfilerStart("/tmp/profile");
#endif

    std::string results_file_path = std::string("/scratch/lfs/pvnick/motif_results/") + get_output_filename();
    MotifFinder engine(SERIES_FILEPATH, TIME_SERIES_LEN, results_file_path);
    engine.run(start_pos, end_pos, K);

#ifdef USE_PROFILER
    ProfilerStop();
#endif
    return 0;
}
