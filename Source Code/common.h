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

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define MIN_RANGE 5
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"

namespace mpi = boost::mpi;
namespace interprocess = boost::interprocess;

extern bool use_mpi;
extern unsigned int query_start_pos();
extern unsigned int search_increment();
extern std::string get_output_filename();
extern std::ostream& msg(std::string str);
extern std::ostream& msgl(std::string str);
extern double lg(double x);


#endif // _COMMON_H_
