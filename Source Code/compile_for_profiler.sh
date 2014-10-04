#!/bin/bash
g++ -I/usr/local/include/gperftools/ -Wall -c -std=c++11 main.cpp
g++ -g -o main main.o -L/usr/local/lib/ -lprofiler 2>&1 | head
