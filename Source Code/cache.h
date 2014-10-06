//todo: I think this is possible to do without shared memory!
//
#ifndef _CACHE_H_
#define _CACHE_H_

#include "common.h"
#include "lemire_envelope.h"
class CacheEntry;
class NonsharedCache;
class SharedCache;

class Cache {
public:
    virtual const CacheEntry& operator[](size_t position) const = 0;
    virtual ~Cache() {};
};

#ifdef USE_MPI
    typedef SharedCache cache_type;
#else
    typedef NonsharedCache cache_type;
#endif

class CacheEntry {
public:
    unsigned int time_series_pos;
    double series_normalized[QUERY_LEN];
    double mean;
    double stddev;
    unsigned int fragment_length;
    LemireEnvelope lemire_envelope;
    CacheEntry(double* time_series, unsigned int position):
        time_series_pos(position),
        mean(0),
        stddev(0)
    {
//todo: double-check the following logic
        double ex = 0, ex2 = 0;
        for (fragment_length = 0; fragment_length + position != TIME_SERIES_LEN && fragment_length != QUERY_LEN; ++fragment_length) {
            double d = time_series[fragment_length + position];
            ex += d;
            ex2 += d*d;
        }
        mean = ex/fragment_length;
        stddev = ex2/fragment_length;
        stddev = sqrt(stddev-mean*mean);
        for(unsigned int i = 0; i != fragment_length; i++) {
             series_normalized[i] = (time_series[i + position] - mean)/stddev;
        }
        lemire_envelope = LemireEnvelope(time_series + position, WARPING_r);
    }
    ~CacheEntry() {}
};

//todo: consider polymorphism here
class NonsharedCache: public Cache {
private:
    typedef std::unique_ptr<CacheEntry[], void(*)(CacheEntry*)> cache_ptr;
    cache_ptr cache;
    static void deallocate_cache(CacheEntry* cache_entries) {
        std::allocator<CacheEntry> cache_alloc;
        for (int i = 0; i != TIME_SERIES_LEN; ++i)
            cache_alloc.destroy(cache_entries + i);
        cache_alloc.deallocate(cache_entries, TIME_SERIES_LEN);
    }
    static bool initialized;
public:
    NonsharedCache() = delete;

    NonsharedCache(double* time_series): cache(nullptr, &NonsharedCache::deallocate_cache) {
        if (initialized) throw std::logic_error("Cache must be initialized only once");
        std::cerr << "Nonshared cache: Allocating " << TIME_SERIES_LEN * sizeof(CacheEntry) << " bytes of memory" << std::endl;
        std::allocator<CacheEntry> cache_alloc;
        std::cerr << "Caching reusable data" << std::endl;
        CacheEntry* tmp_cache = cache_alloc.allocate(TIME_SERIES_LEN);
        for (int i = 0; i != TIME_SERIES_LEN; ++i) {
            cache_alloc.construct(tmp_cache + i, time_series, i);
            if ((i % 100000) == 0) {
                std::cerr << i << "/" << TIME_SERIES_LEN << " (" << ((float)i / TIME_SERIES_LEN * 100) << "% complete)" << std::endl;
            }
        }
        cache.reset(tmp_cache);
        initialized = true;
    }

    const CacheEntry& operator[](size_t position) const {
        return cache[position];
    }
};
bool NonsharedCache::initialized = false;

class SharedCache: public Cache {
    boost::interprocess::mapped_region* cache;
    double* time_series;

    typedef unsigned int rank_id;
    typedef unsigned int rank_id_within_hostname;
    typedef std::map<rank_id, rank_id_within_hostname> hostname_proc_layout;
    typedef std::pair<std::string, rank_id> hostname_rank_pair;

    const rank_id root_rank = 0;
    std::map<std::string, hostname_proc_layout> hostname_proc_layouts;
    bool hostname_proc_layouts_constructed;
    size_t hostname_proc_count;
    std::string global_shm_tag;
    static bool initialized;

    rank_id rank() const {
        mpi::communicator world;
        return world.rank();
    }

    size_t num_procs() const {
        mpi::communicator world;
        return world.size();
    }

    bool is_root() const {
        return rank() == root_rank;
    }

    std::string hostname() const {
        return boost::asio::ip::host_name();
    }

    hostname_proc_layout get_my_hostname_proc_layout() {
        std::string my_hostname = hostname();
        if (hostname_proc_layouts.find(my_hostname) == hostname_proc_layouts.end()) {
            throw std::runtime_error("Hostname not found");
        }
        return hostname_proc_layouts[my_hostname];
    }

    rank_id_within_hostname get_my_intrahostname_rank() {
        hostname_proc_layout my_hostname_layout = get_my_hostname_proc_layout();
        rank_id my_rank = rank();
        if (my_hostname_layout.find(my_rank) == my_hostname_layout.end()) {
            throw std::runtime_error("Rank not found (corrupted hostname layout)");
        }
        return my_hostname_layout[my_rank];
    }

    size_t num_hostname_procs() const {
        return hostname_proc_count;
    }

    void synchronize_random_seed() {
        mpi::communicator world;
        if (is_root()) {
            msgl("Broadcasting random seed");
            time_t seed = time(NULL);
            mpi::broadcast(world, seed, root_rank);
            srand(time(NULL));
        } else {
            time_t seed;
            mpi::broadcast(world, seed, root_rank);
            srand(time(NULL));
        }
        global_shm_tag = std::to_string((unsigned long long)rand());
    }

    void construct_hostname_proc_layouts() {
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
            hostname_proc_count = hostname_proc_layouts[hostname()].size();
            hostname_proc_layouts_constructed = true;
        }
    }

    void sync_hostname_ring() {
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
        mpi::communicator world;
        int token;
        if (my_globalrank_hnrank == my_hostname_layout.begin()) {
            //i am the first proc in the ring
            auto tmp = my_hostname_layout.end();
            prev_proc = (--tmp)->first;
            tmp = my_globalrank_hnrank;
            next_proc = (++tmp)->first;
            token = 1;
            world.send(next_proc, next_proc, token);
            world.recv(prev_proc, rank(), token);
        } else if (my_globalrank_hnrank == --my_hostname_layout.end()) {
            //i am the last proc in the ring
            auto tmp = my_globalrank_hnrank;
            prev_proc = (--tmp)->first;
            tmp = my_hostname_layout.begin();
            next_proc = tmp->first;
            world.recv(prev_proc, rank(), token);
            world.send(next_proc, next_proc, token);
        } else {
            //i am somewhere in the middle of the ring
            auto tmp = my_globalrank_hnrank;
            prev_proc = (--tmp)->first;
            tmp = my_globalrank_hnrank;
            next_proc = (++tmp)->first;
            world.recv(prev_proc, rank(), token);
            world.send(next_proc, next_proc, token);
        }
    }

    std::string shared_fragment_name(size_t fragment_id) const {
        std::ostringstream name;
        name << "Motif_finder_cache_fragment_" << global_shm_tag << "_" << fragment_id;
        msg("looking for ") << name.str() << std::endl;
        return name.str();
    }

    size_t fragment_id_from_timeseries_position(size_t position) const {
        size_t bucket_size = TIME_SERIES_LEN / num_hostname_procs();
        if (position / bucket_size >= num_hostname_procs())
            return num_hostname_procs() - 1;
        return position / bucket_size;
    }

    size_t fragment_position_from_timeseries_position(size_t position) const {
        size_t fragment_id = fragment_id_from_timeseries_position(position);
        size_t bucket_size = TIME_SERIES_LEN / num_hostname_procs();
        return position - fragment_id * bucket_size; 
    }

    void initialize_cache_fragments() {
        //treat all procs in the hostname as existing in a ring
        //a proc is chosen to initialize its shared cache fragment. that proc tells the next proc to
        //initialize its fragment, etc, until the last proc receives the signal and tells the first proc.
        //the first proc then tells its next proc that all cache fragments are initialized, which passes
        //that to the next proc, etc. around the ring.

        size_t my_fragment_id = get_my_intrahostname_rank();
        size_t num_fragments = num_hostname_procs();
        std::string fragment_name = shared_fragment_name(my_fragment_id);

        msg("Allocating shared memory \"") << fragment_name << "\"" << std::endl;
        boost::interprocess::shared_memory_object shm_fragment (
            boost::interprocess::create_only,
            fragment_name.c_str(),
            boost::interprocess::read_write);

        size_t bucket_size = TIME_SERIES_LEN / num_hostname_procs();
        size_t num_entries = TIME_SERIES_LEN - my_fragment_id * bucket_size; //last bucket may hold a few more entries
        shm_fragment.truncate(num_entries * sizeof(CacheEntry));

        std::allocator<boost::interprocess::mapped_region> region_alloc;
        region_alloc.construct(cache + my_fragment_id, shm_fragment, boost::interprocess::read_write);
        CacheEntry* fragment_cache_entries = static_cast<CacheEntry*>(cache[my_fragment_id].get_address());
        /*boost::interprocess::mapped_region region(shm_fragment, boost::interprocess::read_write);
        CacheEntry* fragment_cache_entries = static_cast<CacheEntry*>(region.get_address());*/
        msgl("Initializing cache fragment");
        std::allocator<CacheEntry> cache_alloc;
size_t entries = 0;
        for (size_t i = 0; i != TIME_SERIES_LEN; ++i) {
            if (fragment_id_from_timeseries_position(i) == my_fragment_id) {
                size_t pos = fragment_position_from_timeseries_position(i);
                cache_alloc.construct(fragment_cache_entries + pos, time_series, i);
++entries;
            }
/*            if ((i % 100000) == 0) {
                msg("Fragment ") << i << "/" << TIME_SERIES_LEN << " (" << ((float)i / TIME_SERIES_LEN * 100) << "% complete)" << std::endl;
            }
*/
        }

        msgl("Waiting for hostname ring sync");
        sync_hostname_ring();
        msgl("Ring synced");

        msgl("Fetching all cache fragments for fast reference");
        for (size_t i = 0; i != num_fragments; ++i) {
            if (i != my_fragment_id) { 
                boost::interprocess::shared_memory_object shm_fragment (
                    boost::interprocess::open_only,
                    shared_fragment_name(i).c_str(),
                    boost::interprocess::read_only);
                region_alloc.construct(cache + i,
                    shm_fragment,
                    boost::interprocess::read_only);
            }
        }
        msgl("Shared cache successfully constructed");
    }

public:
    SharedCache() = delete;
    SharedCache(double* series):
        time_series(series),
        hostname_proc_layouts_constructed(false)
    {
        if (initialized) throw std::logic_error("Cache must be initialized only once");
        construct_hostname_proc_layouts();
        synchronize_random_seed();
        std::allocator<boost::interprocess::mapped_region> region_alloc;
        cache = region_alloc.allocate(num_hostname_procs());
        initialize_cache_fragments();
        initialized = true;
    }

    ~SharedCache() {
        std::string my_fragment_name = shared_fragment_name(get_my_intrahostname_rank());
        boost::interprocess::shared_memory_object::remove(my_fragment_name.c_str());
        std::allocator<boost::interprocess::mapped_region> region_alloc;
        for (size_t i = 0; i != num_hostname_procs(); ++i)
            region_alloc.destroy(cache + i);
        region_alloc.deallocate(cache, num_hostname_procs());
    }

    const CacheEntry& operator[](size_t position) const {
        size_t fragment_id = fragment_id_from_timeseries_position(position);
        ptrdiff_t fragment_position = fragment_position_from_timeseries_position(position);
        CacheEntry* entries = static_cast<CacheEntry*>(cache[fragment_id].get_address());
        return entries[fragment_position];
/*
        boost::interprocess::shared_memory_object shm_fragment (
            boost::interprocess::open_only,
            shared_fragment_name(fragment_id).c_str(),
            boost::interprocess::read_only);
        boost::interprocess::mapped_region region(
            shm_fragment,
            boost::interprocess::read_only);
        CacheEntry* entry = static_cast<CacheEntry*>(region.get_address()) + fragment_position;
        return *entry; */
    }

};

bool SharedCache::initialized = false;

#endif // _CACHE_H_
