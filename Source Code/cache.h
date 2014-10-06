#ifndef _CACHE_H_
#define _CACHE_H_
#include "common.h"
#include "lemire_envelope.h"

class CacheEntry;
class NonsharedCache;
class SharedCache;

#ifdef USE_MPI
    typedef SharedCache cache_type;
#else
    typedef NonsharedCache cache_type;
#endif

class CacheEntry {
public:
    unsigned int time_series_pos;
    double series_normalized[QUERY_LEN];
    double* series_window;
    double mean;
    double stddev;
    unsigned int fragment_length;
    LemireEnvelope lemire_envelope;
    CacheEntry(double* time_series, unsigned int position):
        time_series_pos(position),
        series_window(time_series + time_series_pos),
        mean(0),
        stddev(0)
    {
        double ex = 0, ex2 = 0;
        for (fragment_length = 0; fragment_length + position != TIME_SERIES_LEN && fragment_length != QUERY_LEN; ++fragment_length) {
            double d = series_window[fragment_length];
            series_normalized[fragment_length] = d;
            ex += d;
            ex2 += d*d;
        }
        mean = ex/fragment_length;
        stddev = ex2/fragment_length;
        stddev = sqrt(stddev-mean*mean);
        for(unsigned int i = 0; i != fragment_length; i++)
             series_normalized[i] = (series_normalized[i] - mean)/stddev;

        lemire_envelope = LemireEnvelope(time_series + position, WARPING_r);
    }
    ~CacheEntry() {}
};

//todo: consider polymorphism here
class NonsharedCache {
private:
    typedef std::unique_ptr<CacheEntry[], void(*)(CacheEntry*)> cache_ptr;
    static cache_ptr cache;
    static void deallocate_cache(CacheEntry* cache_entries) {
        std::allocator<CacheEntry> cache_alloc;
        for (int i = 0; i != TIME_SERIES_LEN; ++i)
            cache_alloc.destroy(cache_entries + i);
        cache_alloc.deallocate(cache_entries, TIME_SERIES_LEN);
    }
    static bool initialized;
public:
    static void destroy() {}

    static void init(double* series) {
        if ( ! initialized) {
            std::cerr << "Nonshared cache: Allocating " << TIME_SERIES_LEN * sizeof(CacheEntry) << " bytes of memory" << std::endl;
            std::allocator<CacheEntry> cache_alloc;
            std::cerr << "Caching reusable data" << std::endl;
            CacheEntry* tmp_cache = cache_alloc.allocate(TIME_SERIES_LEN);
            for (int i = 0; i != TIME_SERIES_LEN; ++i) {
                cache_alloc.construct(tmp_cache + i, series, i);
                if ((i % 100000) == 0) {
                    std::cerr << i << "/" << TIME_SERIES_LEN << " (" << ((float)i / TIME_SERIES_LEN * 100) << "% complete)" << std::endl;
                }
            }
            cache.reset(tmp_cache);
            initialized = true;
        }
    }

    NonsharedCache() {
        if ( ! initialized) throw std::logic_error("Cache must be initialized before use");
    }

    const CacheEntry& operator[](size_t position) {
        return cache[position];
    }
};

class SharedCache {
    static std::vector<boost::interprocess::mapped_region*> cache;
    static double* time_series;

    typedef unsigned int rank_id;
    typedef unsigned int rank_id_within_hostname;
    typedef std::map<rank_id, rank_id_within_hostname> hostname_proc_layout;
    typedef std::pair<std::string, rank_id> hostname_rank_pair;

    const static rank_id root_rank = 0;
    static std::map<std::string, hostname_proc_layout> hostname_proc_layouts;
    static bool hostname_proc_layouts_constructed;
    static size_t hostname_proc_count;
    static bool initialized;

    static rank_id rank() {
        mpi::communicator world;
        return world.rank();
    }

    static size_t num_procs() {
        mpi::communicator world;
        return world.size();
    }

    static bool is_root() {
        return rank() == root_rank;
    }

    static std::string hostname() {
        return boost::asio::ip::host_name();
    }

    static hostname_proc_layout get_my_hostname_proc_layout() {
        std::string my_hostname = hostname();
        if (hostname_proc_layouts.find(my_hostname) == hostname_proc_layouts.end()) {
            throw std::runtime_error("Hostname not found");
        }
        return hostname_proc_layouts[my_hostname];
    }

    static rank_id_within_hostname get_my_intrahostname_rank() {
        hostname_proc_layout my_hostname_layout = get_my_hostname_proc_layout();
        rank_id my_rank = rank();
        if (my_hostname_layout.find(my_rank) == my_hostname_layout.end()) {
            throw std::runtime_error("Rank not found (corrupted hostname layout)");
        }
        return my_hostname_layout[my_rank];
    }

    static size_t num_hostname_procs() {
        return hostname_proc_count;
    }

    static void construct_hostname_proc_layouts() {
        mpi::communicator world;
        if ( ! hostname_proc_layouts_constructed) {
            msgl("Constructing hostname proc layouts");
            hostname_rank_pair my_hostname_rank(hostname(), rank());
            if (is_root()) {
                std::vector<hostname_rank_pair> all_hostname_ranks;
                //root gathers hostname-rank pairs from all procs, assigned intrahostname ranks,
                //then broadcasts the info out to all procs
                mpi::gather(world, my_hostname_rank, all_hostname_ranks, root_rank);
                std::map<std::string, rank_id_within_hostname> hostname_proc_counts;
                for (auto hostname_rank: all_hostname_ranks) {
                    std::string hn = hostname_rank.first;
                    rank_id r = hostname_rank.second;
                    hostname_proc_layouts[hn][r] = hostname_proc_counts[hn]++;
                }
                msgl("I am the root node and I have assigned the following intrahostname ranks:");
                for (auto layout: hostname_proc_layouts) {
                    msg("Hostname ") << layout.first << ":" << std::endl;
                    for (auto assignment: layout.second) {
                        msg("    Proc ") << assignment.first << ": rank" << assignment.second << std::endl;
                    }
                }
                msgl("Broadcasting roster");
                mpi::broadcast(world, hostname_proc_layouts, root_rank);
            } else {
                //this proc is not the root, so send the hostname-rank pair to the root, then wait
                //for the root to broadcast the rank assignments
                mpi::gather(world, my_hostname_rank, root_rank);
                mpi::broadcast(world, hostname_proc_layouts, root_rank);
                msg("I have been assigned rank ") << get_my_intrahostname_rank() << " among host " << hostname() << std::endl;
            }
            hostname_proc_count = hostname_proc_layouts.size();
            hostname_proc_layouts_constructed = true;
        }
    }

    static void sync_hostname_ring() {
        //treat all procs in the hostname as existing in a ring
        //a proc is chosen as the first proc. that proc sends a signal to its next proc, etc,
        //all around the ring until the last proc receives the signal and tells the first proc
        //this allows everyone on the host to sync to a common point
        hostname_proc_layout my_hostname_layout = get_my_hostname_proc_layout();
        auto globalrank_hnrank = my_hostname_layout.begin();
        for (; globalrank_hnrank->first != rank(); ++globalrank_hnrank);
        auto my_globalrank_hnrank = globalrank_hnrank;
        rank_id prev_proc;
        rank_id next_proc;
        if (my_globalrank_hnrank == my_hostname_layout.begin()) {
            //i am the first proc in the ring
            auto tmp = my_hostname_layout.end();
            prev_proc = (--tmp)->first;
            tmp = my_globalrank_hnrank;
            next_proc = (++tmp)->first;
        } else if (my_globalrank_hnrank == --my_hostname_layout.end()) {
            //i am the last proc in the ring
            auto tmp = my_globalrank_hnrank;
            prev_proc = (--tmp)->first;
            tmp = my_hostname_layout.begin();
            next_proc = tmp->first;
        } else {
            //i am somewhere in the middle of the ring
            auto tmp = my_globalrank_hnrank;
            prev_proc = (--tmp)->first;
            tmp = my_globalrank_hnrank;
            next_proc = (++tmp)->first;
        }

        mpi::communicator world;
        mpi::request reqs[2];
        bool send_val = true, recv_val;
        reqs[0] = world.isend(next_proc, rank(), send_val);
        reqs[1] = world.irecv(prev_proc, rank(), recv_val);
        mpi::wait_all(reqs, reqs + 2);
    }

    static const char* shared_fragment_name(size_t fragment_id) {
        std::ostringstream name;
        name << "Motif finder cache fragment " << fragment_id;
        return name.str().c_str();
    }

    static size_t fragment_id_from_timeseries_position(size_t position) {
        return position / num_hostname_procs();
    }

    static size_t fragment_position_from_timeseries_position(size_t position) {
        return position % num_hostname_procs();
    }

    static void initialize_cache_fragments() {
        //treat all procs in the hostname as existing in a ring
        //a proc is chosen to initialize its shared cache fragment. that proc tells the next proc to
        //initialize its fragment, etc, until the last proc receives the signal and tells the first proc.
        //the first proc then tells its next proc that all cache fragments are initialized, which passes
        //that to the next proc, etc. around the ring.

        size_t my_fragment_id = get_my_intrahostname_rank();
        size_t num_fragments = num_hostname_procs();
        const char* fragment_name = shared_fragment_name(my_fragment_id);

        msg("Creating shared memory object \"") << fragment_name << "\"" << std::endl;
        boost::interprocess::shared_memory_object::remove(fragment_name);
        boost::interprocess::shared_memory_object my_fragment_shm(
            boost::interprocess::create_only,
            fragment_name,
            boost::interprocess::read_write);

        msgl("Allocating shared cache fragment memory");
        my_fragment_shm.truncate(TIME_SERIES_LEN * sizeof(CacheEntry) / num_fragments + 1);
        cache[my_fragment_id] = new interprocess::mapped_region(my_fragment_shm, boost::interprocess::read_write);
        CacheEntry* fragment_cache_entries = static_cast<CacheEntry*>(cache[my_fragment_id]->get_address());

        msgl("Initializing cache fragment");
        std::allocator<CacheEntry> cache_alloc;
        for (size_t i = 0; i != TIME_SERIES_LEN; ++i) {
            if (fragment_id_from_timeseries_position(i) == my_fragment_id) {
                size_t pos = fragment_position_from_timeseries_position(i);
                cache_alloc.construct(fragment_cache_entries + pos, time_series, i);
            }
            if ((i % 100000) == 0) {
                std::cerr << i << "/" << TIME_SERIES_LEN << " (" << ((float)i / TIME_SERIES_LEN * 100) << "% complete)" << std::endl;
            }
        }

        msgl("Waiting for hostname ring sync");
        sync_hostname_ring();
        msgl("Ring synced");

        msgl("Fetching all cache fragments for fast reference");
        for (size_t i = 0; i != num_fragments; ++i) {
            if (i != my_fragment_id) {
                boost::interprocess::shared_memory_object fragment_shm(
                    boost::interprocess::open_only,
                    shared_fragment_name(i),
                    boost::interprocess::read_only);

                cache[my_fragment_id] = new interprocess::mapped_region(fragment_shm, boost::interprocess::read_only);
            }
        }

        msgl("Shared cache successfully constructed");
    }
public:
    static void init(double* series) {
        if ( ! initialized) {
            time_series = series;
            construct_hostname_proc_layouts();
            cache = std::vector<boost::interprocess::mapped_region*>(num_hostname_procs());
            initialize_cache_fragments();
            initialized = true;
        }
    }

    static void destroy() {
        const char* my_fragment_name = shared_fragment_name(get_my_intrahostname_rank());
        boost::interprocess::shared_memory_object::remove(my_fragment_name);
    }

    SharedCache() {
        if ( ! initialized) throw std::logic_error("Shared cache must be initialized before use");
    }

    const CacheEntry& operator[](size_t position) {
        size_t fragment_id = fragment_id_from_timeseries_position(position);
        boost::interprocess::mapped_region* fragment = cache[fragment_id];
        CacheEntry* fragment_cache_entries = static_cast<CacheEntry*>(fragment->get_address());

        size_t fragment_position = fragment_position_from_timeseries_position(position);
        return fragment_cache_entries[fragment_position];
    }

};

#endif // _CACHE_H_
