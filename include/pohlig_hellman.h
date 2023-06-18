#include <cstdint>
#include <gmpxx.h>
#include <vector> 

namespace PohligHellman {
    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p, std::vector<std::pair<mpz_class, mpz_class>> &factors);
};