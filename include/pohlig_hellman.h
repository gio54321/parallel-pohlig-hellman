#include <cstdint>
#include <gmpxx.h>
#include <vector> 

namespace PohligHellman {
    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p, std::vector<std::pair<mpz_class, mpz_class>> &factors);
    mpz_class discrete_log_parallel(const mpz_class g, const mpz_class b, const mpz_class p, const std::vector<std::pair<mpz_class, mpz_class>> &factors, int num_workers, int num_workers_bsgs);
};