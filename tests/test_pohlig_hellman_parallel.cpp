#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <sstream>
#include <cassert>
#include <gmpxx.h>

#include "discrete_utils.h"
#include "utimer.h"
#include "pohlig_hellman.h"

void test_pohlig_hellman(int num_workers, int num_bsgs_workers) {
    std::cout << "num_workers: " << num_workers << std::endl;
    std::ifstream myfile;
    myfile.open("input/pohlig_hellman_inputs.txt");

    std::string line;
    if (myfile.is_open()) {
        while(std::getline(myfile, line)){
            std::cout << line << std::endl;
            
            std::getline(myfile, line);
            std::string sg, sb, sp;
            std::istringstream(line) >> sg >> sb >> sp;
            mpz_class g(sg);
            mpz_class b(sb);
            mpz_class p(sp);

            // std::cout << "g: " << g << ", b: " << b << ", p: " << p << std::endl;

            std::getline(myfile, line);
            auto line_stream = std::istringstream(line);

            int size = 0;
            line_stream >> size;

            std::vector<std::pair<mpz_class, mpz_class>> factorization;
            for (int i = 0; i < size; ++i) {
                std::string s_prime, s_e;
                line_stream >> s_prime >> s_e;
                mpz_class prime(s_prime);
                mpz_class e(s_e);
                factorization.push_back({prime, e});
            }

            // print factorization
            // std::cout << "factorization: ";
            // for (auto factor : factorization) {
            //     std::cout << "(" << factor.first << ", " << factor.second << ") ";
            // }
            // std::cout << std::endl;

            long time_taken;
            mpz_class result;
            {
                utimer u("discrete log", &time_taken);
                if (num_workers == 0) {
                    result = PohligHellman::discrete_log(g, b, p, factorization);
                } else {
                    result = PohligHellman::discrete_log_parallel(g, b, p, factorization, num_workers, num_bsgs_workers);
                }
            }

            // std::cout << "result: " << result << std::endl;
            std::cout << time_taken << " us" << std::endl;


            // check that we actually found the discrete log
            assert(powerMod(g, result, p) == b);

            std::getline(myfile, line);
        }
    }
}




int main(int argc, char** argv) {
    if (argc != 3) {
        std::cout << "Usage: ./main <num_workers> <num_bsgs_workers>" << std::endl;
        std::cout << "num_workers: number of workers to use for Pohlig-Hellman" << std::endl;
        std::cout << "num_bsgs_workers: number of workers to use for BSGS" << std::endl;
        std::cout << "total number of workers used: num_workers * num_bsgs_workers" << std::endl;
        std::cout << "if num_workers = 0 then use sequential version" << std::endl;
        return 1;
    } 
    test_pohlig_hellman(atoi(argv[1]), atoi(argv[2]));

    return 0;
}