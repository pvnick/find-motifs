#include "find_motifs.h"
#include <iostream>
#include <string>


/// Main Function
int main(  int argc , char *argv[] )
{
    unsigned int K = 100;
    MotifFinder engine;
    for (int i = 0; i != 100; ++i) {
        MotifFinder::TopKMatches result = engine.single_pass(K, i);
        result.print();
    }
    return 0;
}
