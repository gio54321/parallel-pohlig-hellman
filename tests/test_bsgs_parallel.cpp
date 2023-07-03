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


void test_bsgs(int num_workers) {
    std::ifstream myfile;
    myfile.open("input/bsgs_inputs.txt");

    int i = 0;
    long time_sum = 0;
    long time_max = 0;
    long time_min = 0;
    std::string line;
    if (myfile.is_open()) {
        while(std::getline(myfile, line)) {
            std::string sg, sb, sp;
            std::istringstream(line) >> sg >> sb >> sp;
            mpz_class g(sg);
            mpz_class b(sb);
            mpz_class p(sp);

            // std::cout << "g: " << g << ", b: " << b << ", p: " << p << std::endl;

            long time_taken;
            mpz_class result;
            {
                utimer u("check", &time_taken);
                result = BabyStepGiantStep::discrete_log_parallel(g, b, p, p-1, num_workers);
            }
            // std::cout << "result: " << result << std::endl;

            time_sum += time_taken;
            if (i == 0) {
                time_max = time_taken;
                time_min = time_taken;
            } else {
                if (time_taken > time_max) {
                    time_max = time_taken;
                }
                if (time_taken < time_min) {
                    time_min = time_taken;
                }
            }

            // check that we actually found the discrete log
            assert(powerMod(g, result, p) == b);

            if (i == 4) {
                i = 0;
                std::cout << p << ", " << (double)time_sum/5 << ", " << time_min << ", " << time_max << std::endl;
                time_sum = 0;
            } else {
                i++;
            }
        }
    }
}




int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: ./main <num_workers>" << std::endl;
        return 1;
    } 
    test_bsgs(atoi(argv[1]));
    return 0;
}