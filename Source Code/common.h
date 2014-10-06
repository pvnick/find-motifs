#ifndef _COMMON_H_
#define _COMMON_H_

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

extern bool use_mpi;
extern std::string global_shm_tag;
extern unsigned int query_start_pos();
extern unsigned int search_increment();
extern std::string get_output_filename();
extern std::ostream& msg(std::string str);
extern std::ostream& msgl(std::string str);

#endif // _COMMON_H_
