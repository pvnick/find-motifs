#define USE_MPI
//#define USE_PROFILER

#include "find_motifs.h"
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>

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
        std::cout << "Proc " << world.rank() << ": " << str;
        return std::cout;
    }

    std::string hostname() {
        return boost::asio::ip::host_name();
    }

    bool is_root() {
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

    rank_id get_hostname_leader() {
        elect_hostname_leaders();
        rank_id hostname_leader = hostname_leaders[hostname()];
        msg("Leader among ") << hostname() << " is " << hostname_leader << std::endl;
        return hostname_leader;
    }

    bool is_hostname_leader() {
        mpi::communicator world;
        return world.rank() == get_hostname_leader();
    }

    MotifFinder get_engine() {
        //engine cache takes up a few gigs of memory. therefore, only initialize it once per machine,
        //then share the memory among all other procs on the same machine
        mpi::communicator world;
        rank_id my_rank = world.rank();
        if (is_hostname_leader()) {
            //current process is the leader on this machine. tell the motif finder engine to allocate the cache, then share its memory
            msg("Initializing cache") << std::endl;
            MotifFinder engine(true, true);
            bool initialized = true;
            mpi::broadcast(world, initialized, hostname_leader);
        } else {
            //a different process is the leader on this machine. wait for it to allocate the cache then tell
            //the motif finder engine to use its shared memory
            msg("Waiting for cache to initialize") << std::endl;
            bool initialized;
            mpi::broadcast(world, initialized, hostname_leader);
            msg("Found initialized cache") << std::endl;
            MotifFinder engine(true, false);
        }

        return engine;
    }

    size_t query_start_position_by_params(size_t series_length, rank_id proc_rank, unsigned int num_procs) {
        proc_rank = num_procs - proc_rank; //make assignment easier to conceptualize
        size_t search_space_length = floor((float)proc_rank / num_procs * series_length * (series_length + 1.0) / 2.0);
        size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
        return position;
    }

    size_t my_query_start_pos() {
        mpi::communicator world;
        unsigned int num_procs = world.size();
        rank_id my_rank = world.rank();
        return query_start_position_by_params(TIME_SERIES_LEN, my_rank, num_procs);
    }

    size_t my_query_end_pos() {
        mpi::communicator world;
        unsigned int num_procs = world.size();
        rank_id my_rank = world.rank();
        return query_start_position_by_params(TIME_SERIES_LEN, my_rank + 1, num_procs);
    }
#else
    size_t my_query_start_pos() {
        return 0;
    }

    size_t my_query_end_pos() {
        return TIME_SERIES_LEN;
    }

    MotifFinder get_engine() {
        MotifFinder engine(false, false);
        return engine;
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
        std::cout << "Start: " << start << ", End: " << fin << ", Queries: " << fin - start << ", Candidates: ";
        std::cout << (TIME_SERIES_LEN - start) * (TIME_SERIES_LEN - start + 1) / 2 - (TIME_SERIES_LEN - fin) * (TIME_SERIES_LEN - fin + 1) / 2 << std::endl;
    }
    */
    unsigned int K = 100;
    size_t start_pos = my_query_start_pos();
    size_t end_pos = my_query_end_pos();
    MotifFinder engine = get_engine();
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
