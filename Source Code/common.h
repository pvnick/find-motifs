#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <sstream>
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

typedef double dist_type;

extern bool use_mpi;
extern unsigned int query_start_pos();
extern unsigned int search_increment();
extern std::string get_output_filename();
extern std::ostream& msg(std::string str);
extern std::ostream& msgl(std::string str);

double lg(double x) {
    return log(x) / log(2);
}

#endif // _COMMON_H_
