#!/bin/bash

module load boost/1.53.0
module load intel/2013_sp1.3.174
module load openmpi/1.8.1

INPUT_FILES="main.cpp"
OUTPUT_FILE="main_mpi"

/apps/mpi/intel/2013/sp1.3.174/openmpi/1.8.1/bin/mpic++ -g -std=c++11 -I/apps/boost/1.53.0/include/ -I/apps/mpi/intel/2013/sp1.3.174/openmpi/1.8.1/include/ $INPUT_FILES -o $OUTPUT_FILE -L/apps/boost/1.53.0/stage/lib/ -lboost_mpi -lboost_serialization -lboost_system -lrt
