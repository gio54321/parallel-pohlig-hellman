#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <sstream>
#include <cassert>
#include <gmpxx.h>

#include "discrete_utils.h"
#include "bsgs.h"
#include "utimer.h"


int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Usage: ./main <input_filename> <num_workers> <load_factor>" << std::endl;
        std::cout << "if num_workers == 0, then use the sequential version" << std::endl;
        std::cout << "if num_workers == -1, then use the \"old\" sequential version that uses std::unordered_map" << std::endl;
        return 1;
    } 

    char *input_filename = argv[1];
    int num_workers = atoi(argv[2]);
    float load_factor = atof(argv[3]);
    
    // read from the input file g, b, p
    std::ifstream myfile;
    myfile.open(input_filename);
    if (!myfile.is_open()) {
        std::cout << "Could not open file: " << input_filename << std::endl;
        return 1;
    }
    std::string line;
    getline(myfile, line);
    std::string sg, sb, sp;
    std::istringstream(line) >> sg >> sb >> sp;
    mpz_class g(sg);
    mpz_class b(sb);
    mpz_class p(sp);

    // run the discrete log algorithm
    long time_taken;
    mpz_class result;
    {
        utimer u("check", &time_taken);
        if (num_workers == 0) {
            result = BabyStepGiantStep::discrete_log(g, b, p, p-1);
        } else if (num_workers == -1) {
            result = BabyStepGiantStep::discrete_log_old(g, b, p, p-1);
        } else {
            result = BabyStepGiantStep::discrete_log_parallel(g, b, p, p-1, num_workers, load_factor);
        }
    }

    // assure that we actually found the discrete log
    assert(powerMod(g, result, p) == b);

    std::cout << "Found result in: " << time_taken << " microseconds" << std::endl;
    return 0;
}