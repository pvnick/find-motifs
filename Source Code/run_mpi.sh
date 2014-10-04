#!/bin/bash

module load boost/1.53.0
module load intel/2013_sp1.3.174
module load openmpi/1.8.1

mpirun -np 8 main_mpi
