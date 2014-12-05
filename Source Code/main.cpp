//config
#define USE_MPI
//#define USE_PROFILER

#include "common.h"
#include "find_motifs.h"
#include "ucr_dtw.h"
#include "cli_options.h"
#include "postprocessor.h"
#include <boost/regex.hpp>

#ifdef USE_MPI
    bool use_mpi = true;
#else
    bool use_mpi = false;
#endif

namespace program_options = boost::program_options;

double lg(double x) {
    return log(x) / log(2);
}

unsigned int get_query_start_pos() {
    if (use_mpi) {
        mpi::communicator world;
        return world.rank();
    } else {
        return 0;
    }
}

unsigned int get_query_increment() {
    if (use_mpi) {
        mpi::communicator world;
        return world.size();
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
    #include <profiler.h>
#endif


/// Main Function
int main(int argc, char *argv[])
{
#ifdef USE_MPI
    mpi::environment env(argc, argv);
#endif
    CLIOptions::init(argc, argv);
    CLIOptions opts = CLIOptions::get_instance();

    std::string command = opts["command"].as<std::string>();
    std::string results_file_path = std::string("/scratch/lfs/pvnick/motif_results/") + get_output_filename();
    std::string series_filepath = "/home/pvnick/oximetry/data/oximetry.txt";
    unsigned int K = opts["K"].as<unsigned int>();
    MotifFinder engine(series_filepath, TIME_SERIES_LEN, K, results_file_path, '\t');

    if (command == "postprocess") {
        SubsequenceLookup& subsequences = engine.get_subsequence_lookup();
        if (opts.count("input-files")) {
            std::vector<std::string> input_files = opts["input-files"].as<std::vector<std::string>>();
            PostProcess::PostProcessor post_processor(input_files, &subsequences, QUERY_LEN);
            post_processor.run();
        } else {
            std::cerr << "--input-files required when using postprocess" << std::endl;
            return 1;
        }
    } else if (command == "find-motifs") {
        size_t query_start_pos = get_query_start_pos();
        size_t query_increment = get_query_increment();
#ifdef USE_PROFILER
        ProfilerStart("/tmp/profile");
#endif

        //engine.run(query_start_pos, query_increment);
        std::vector<double> const& time_series = engine.get_timeseries();
        std::cout << time_series[2] << std::endl;
        SubsequenceLookup subsequences = engine.get_subsequence_lookup();
        Subsequence const& subsequence = subsequences[5];
        std::cout << subsequence.time_series_pos << std::endl;

#ifdef USE_PROFILER
        ProfilerStop();
#endif
    } else {
        std::cerr << "Unknown command " << command << std::endl;
        return 1;
    }

    return 0;
}
