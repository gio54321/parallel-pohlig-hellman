#ifndef BSGS_H
#define BSGS_H

#include <cstdint>
#include <gmpxx.h>

namespace BabyStepGiantStep {
    const float default_load_factor = 4.5;

    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p);
    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p, mpz_class order);
    mpz_class discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers);
    mpz_class discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers, float load_factor);
};

#endif