#define USE_MPI
//#define USE_PROFILER

#include "find_motifs.h"
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <sstream>
#include <fstream>

#ifdef USE_PROFILER
    #include "profiler.h"
#endif
#ifdef USE_MPI
    #include <boost/mpi.hpp>
    #include <boost/serialization/string.hpp>
    #include <boost/serialization/map.hpp>
    #include <boost/asio.hpp>

    namespace mpi = boost::mpi;
    namespace interprocess = boost::interprocess;

    typedef unsigned int rank_id;
    typedef std::pair<std::string, rank_id> hostname_rank_pair;
    static rank_id root_rank = 0;
    static std::map<std::string, rank_id> hostname_leaders;
    static bool hostname_leaders_elected = false;

    std::ostream& msg(std::string str) {
        mpi::communicator world;
        std::cerr << "Proc " << world.rank() << ": " << str;
        return std::cerr;
    }

    const std::string hostname() {
        return boost::asio::ip::host_name();
    }

    const bool is_root() {
        mpi::communicator world;
        return world.rank() == root_rank;
    }

    void elect_hostname_leaders() {
        mpi::communicator world;
        if ( ! hostname_leaders_elected) {
            hostname_rank_pair my_hostname_rank(hostname(), world.rank());
            if (is_root()) {
                std::vector<hostname_rank_pair> all_hostname_ranks;
                //root gathers hostname-rank pairs from all procs, elects the hostname leaders,
                //then broadcasts the leader info out to all procs
                mpi::gather(world, my_hostname_rank, all_hostname_ranks, root_rank);
                for (auto hostname_rank: all_hostname_ranks) {
                    std::string hn = hostname_rank.first;
                    rank_id r = hostname_rank.second;
                    if (hostname_leaders.find(hn) == hostname_leaders.end())
                        hostname_leaders[hn] = r;
                    else
                        hostname_leaders[hn] = std::min(hostname_leaders[hn], r); //leader rank is minimum rank for that hostname
                }
                mpi::broadcast(world, hostname_leaders, root_rank);
            } else {
                //this proc is not the root, so send the hostname-rank pair to the root, then wait
                //for the root to broadcast the leaders
                mpi::gather(world, my_hostname_rank, root_rank);
                mpi::broadcast(world, hostname_leaders, root_rank);
            }
            hostname_leaders_elected = true;
        }
    }

    const rank_id get_hostname_leader() {
        elect_hostname_leaders();
        rank_id hostname_leader = hostname_leaders[hostname()];
        msg("Leader among ") << hostname() << " is " << hostname_leader << std::endl;
        return hostname_leader;
    }

    const bool is_hostname_leader() {
        mpi::communicator world;
        return world.rank() == get_hostname_leader();
    }

    MotifFinder* get_engine() {
        //engine cache takes up a few gigs of memory. therefore, only initialize it once per machine,
        //then share the memory among all other procs on the same machine
        MotifFinder* engine;
        mpi::communicator world;
        rank_id my_rank = world.rank();
        if (is_hostname_leader()) {
            //current process is the leader on this machine. tell the motif finder engine to allocate the cache, then share its memory
            msg("Initializing cache") << std::endl;
            engine = new MotifFinder(true, true);
            bool initialized = true;
            mpi::broadcast(world, initialized, get_hostname_leader());
        } else {
            //a different process is the leader on this machine. wait for it to allocate the cache then tell
            //the motif finder engine to use its shared memory
            msg("Waiting for cache to initialize") << std::endl;
            bool initialized;
            mpi::broadcast(world, initialized, get_hostname_leader());
            msg("Found initialized cache") << std::endl;
            engine = new MotifFinder(true, false);
        }
        return engine;
    }

    const size_t query_start_position_by_params(size_t series_length, rank_id proc_rank, unsigned int num_procs) {
        proc_rank = num_procs - proc_rank; //make assignment easier to conceptualize
        size_t search_space_length = floor((float)proc_rank / num_procs * series_length * (series_length + 1.0) / 2.0);
        size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
        return position;
    }

    const size_t my_query_start_pos() {
        mpi::communicator world;
        unsigned int num_procs = world.size();
        rank_id my_rank = world.rank();
        return query_start_position_by_params(TIME_SERIES_LEN, my_rank, num_procs);
    }

    const size_t my_query_end_pos() {
        mpi::communicator world;
        unsigned int num_procs = world.size();
        rank_id my_rank = world.rank();
        return query_start_position_by_params(TIME_SERIES_LEN, my_rank + 1, num_procs);
    }

    const std::string get_output_filename() {
        std::ostringstream filename;
        filename << "results" << world.rank() << ".out";
        return filename.c_str();
    }
#else
    std::ostream& msg(std::string str) {
        std::cerr << str;
        return std::cerr;
    }

    const size_t my_query_start_pos() {
        return 0;
    }

    const size_t my_query_end_pos() {
        return TIME_SERIES_LEN;
    }

    MotifFinder* get_engine() {
        return new MotifFinder(false, false);
    }

    const std::string get_output_filename() {
        return "results.out";
    }

#endif

/// Main Function
int main(int argc , char *argv[] )
{
#ifdef USE_MPI
    mpi::environment env(argc, argv);
#endif
    /*
    int num_procs = 100;
    for (int i = 0; i != num_procs; ++i) {
        size_t fin = query_start_position(TIME_SERIES_LEN, i + 1, num_procs);
        size_t start = query_start_position(TIME_SERIES_LEN, i, num_procs);
        std::cerr << "Start: " << start << ", End: " << fin << ", Queries: " << fin - start << ", Candidates: ";
        std::cerr << (TIME_SERIES_LEN - start) * (TIME_SERIES_LEN - start + 1) / 2 - (TIME_SERIES_LEN - fin) * (TIME_SERIES_LEN - fin + 1) / 2 << std::endl;
    }
    */
    unsigned int K = 100;
    size_t start_pos = my_query_start_pos();
    size_t end_pos = my_query_end_pos();
    MotifFinder* engine = get_engine();
#ifdef USE_PROFILER
    ProfilerStart("/tmp/profile");
    end_pos = start_pos + 3;
#endif
    std::ofstream out(std::string("/scratch/lfs/pvnick/motif_results/") + get_output_filename());
    for (size_t i = start_pos; i != end_pos
                            && i != TIME_SERIES_LEN - QUERY_LEN /*don't query at the end of the time series*/
                            ; ++i) {
        msg("Querying from ") << i << std::endl;
        MotifFinder::TopKMatches results = engine->single_pass(K, i);
        results >> out;
        msg("Finished ") << i << "/" << (end_pos - start_pos) << " (" << (float)i / (end_pos - start_pos) * 100 << "%) complete" << std::endl;
    }
    out.close();
#ifdef USE_PROFILER
    ProfilerStop();
#endif
    delete engine;
    return 0;
}
