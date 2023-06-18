#include <cstdint>
#include <gmpxx.h>

namespace BabyStepGiantStep {
    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p);
    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p, mpz_class order);
    uint64_t discrete_log_native(uint64_t g, uint64_t b, uint64_t p);
    uint64_t discrete_log_native_asm(uint64_t g, uint64_t b, uint64_t p);
};