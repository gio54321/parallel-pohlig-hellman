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

mpz_class dlog_algo(mpz_class g, mpz_class b, mpz_class p, int algo, long *time_taken) {

    mpz_class result;
    {
        utimer u("discrete log", time_taken);
        BabyStepGiantStep bsgs;
        if (algo == 1) {
            result = bsgs.discrete_log(g, b, p);
        } else if (algo == 2) {
            result = bsgs.discrete_log_2(g, b, p);
        } else if (algo == 3) {
            result = bsgs.discrete_log_native(g.get_ui(), b.get_ui(), p.get_ui());
        } else if (algo == 4) {
            result = bsgs.discrete_log_native_asm(g.get_ui(), b.get_ui(), p.get_ui());
        }
    }
    return result;
}

void test_bsgs(int algo) {
    if (algo == 1) {
        std::cout << "algo: bsgs" << std::endl;
    } else if (algo == 2) {
        std::cout << "algo: bsgs_2" << std::endl;
    } else if (algo == 3) {
        std::cout << "algo: bsgs_native" << std::endl;
    } else if (algo == 4) {
        std::cout << "algo: bsgs_native_asm" << std::endl;
    }

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

            long time_taken;
            mpz_class result = dlog_algo(g, b, p, algo, &time_taken);

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
        std::cout << "Usage: ./main <algo>" << std::endl;
        std::cout << "algo: 1 - bsgs, 2 - bsgs_2, 3 - bsgs_native, 4 - bsgs_native_asm" << std::endl;
        return 1;
    } 
    test_bsgs(atoi(argv[1]));
    return 0;
}