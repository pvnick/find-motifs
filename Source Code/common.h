#ifndef _COMMON_H_
#define _COMMON_H_

#define USE_MPI
//#define USE_PROFILER
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iterator>

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"

#ifdef USE_MPI
    #include <boost/mpi.hpp>
    #include <boost/interprocess/shared_memory_object.hpp>
    #include <boost/interprocess/mapped_region.hpp>
    #include <boost/serialization/string.hpp>
    #include <boost/serialization/map.hpp>
    #include <boost/asio.hpp>
    namespace mpi = boost::mpi;
    namespace interprocess = boost::interprocess;

    unsigned int query_start_pos() {
        mpi::communicator world;
        return world.rank();
    }

    unsigned int search_increment() {
        mpi::communicator world;
        return world.size();
    }

    std::string get_output_filename() {
        mpi::communicator world;
        std::ostringstream filename;
        filename << "results" << world.rank() << ".out";
        return filename.str();
    }

    std::ostream& msg(std::string str) {
        mpi::communicator world;
        std::cerr << "Proc " << world.rank() << ": " << str;
        return std::cerr;
    }

    std::ostream& msgl(std::string str) {
        return msg(str) << std::endl;
    }

#else
    unsigned int query_start_pos() {
        return 0;
    }

    unsigned int search_increment() {
        return 1;
    }

    std::string get_output_filename() {
        return "results.out";
    }

    std::ostream& msg(std::string str) {
        std::cerr << str;
        return std::cerr;
    }

    std::ostream& msgl(std::string str) {
        return msg(str) << std::endl;
    }
#endif //USE_MPI

#endif //_COMMON_H_
