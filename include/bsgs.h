#ifndef BSGS_H
#define BSGS_H

#include <cstdint>
#include <gmpxx.h>

namespace BabyStepGiantStep {
    const float default_load_factor = 4.5;

    mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p, mpz_class order);

    inline mpz_class discrete_log(mpz_class g, mpz_class b, mpz_class p) {
        mpz_class order = p - 1;
        return discrete_log(g, b, p, order);
    }

    mpz_class discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers, float load_factor);

    inline mpz_class discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers) {
        return discrete_log_parallel(g, b, p, order, num_workers, default_load_factor);
    }

    // old sequential implementation, left for reference
    mpz_class discrete_log_old(mpz_class g, mpz_class b, mpz_class p, mpz_class order);
};

#endif