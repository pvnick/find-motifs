#include "common.h"
#include "shared_cache.h"

bool NonsharedCache::initialized = false;
NonsharedCache::cache_ptr NonsharedCache::cache = cache_ptr(nullptr, NonsharedCache::deallocate_cache);

#ifdef USE_MPI
    bool SharedCache::hostname_proc_layouts_constructed = false;
    size_t SharedCache::hostname_proc_count = 0;
    bool SharedCache::initialized = false;
#endif
