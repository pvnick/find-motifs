#include "find_motifs.h"
#include <iostream>
#include <string>


/// Main Function
int main(  int argc , char *argv[] )
{
    std::string time_series_file(argv[1]);
    unsigned int time_series_length(atoi(argv[2]));
    double warping_window(atof(argv[3]));
    unsigned int K = 100;
    MotifFinder engine(time_series_file, time_series_length, 100, warping_window);
    std::cout << (unsigned long long)time_series_length * time_series_length << std::endl;
    for (int i = 0; i != 100; ++i) {
        MotifFinder::TopKMatches result = engine.single_pass(K, i);
        result.print();
    }
    return 0;
}
