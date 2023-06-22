#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cstdint>
#include <gmpxx.h>
#include <vector>
#include <thread>

#include "pohlig_hellman.h"
#include "bsgs.h"
#include "discrete_utils.h"


// compute discrete log for a group of prime power order
// https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
mpz_class discrete_log_prime_power(mpz_class g, mpz_class b, mpz_class mod, mpz_class p, mpz_class e)
{
    //std::cout << "discrete_log_prime_power" << std::endl;
    //std::cout << "g: " << g << ", b: " << b << ", mod: " << mod << ", p: " << p << ", e: " << e << std::endl;

    if (e == 1) {
        return BabyStepGiantStep::discrete_log(g, b, mod, p);
    }
    mpz_class result = 0;

    const mpz_class gamma = powerMod(g, powerMod(p, e-1, mod), mod);

    for (mpz_class k = 0; k < e; ++k) {
        const mpz_class exponent = powerMod(p, e-1-k, mod);
        const mpz_class h_k = powerMod((b * powerMod(modInversePrime(g, mod), result, mod)) % mod, exponent, mod);
        const mpz_class dk = BabyStepGiantStep::discrete_log(gamma, h_k, mod, p);

        result += (dk * powerMod(p, k, mod)) % mod;
    }
    return result;
}

// compute discrete log with Pohlig-Hellman algorithm
// https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
mpz_class PohligHellman::discrete_log(mpz_class g, mpz_class b, mpz_class p, std::vector<std::pair<mpz_class, mpz_class>> &factors)
{
    const mpz_class order = p - 1;

    std::vector<mpz_class> x;
    x.reserve(factors.size());

    std::vector<mpz_class> h_i;
    h_i.reserve(factors.size());

    // iterate over the factors
    for (auto factor : factors) {
        const mpz_class prime = factor.first;
        const mpz_class e = factor.second;

        const mpz_class p_to_e = powerMod(prime, e, p);
        const mpz_class subgroup_order = order / p_to_e;

        const mpz_class subgroup_gen = powerMod(g, subgroup_order, p);
        const mpz_class h = powerMod(b, subgroup_order, p);

        const mpz_class x_i = discrete_log_prime_power(subgroup_gen, h, p, prime, e);

        x.push_back(x_i);
        h_i.push_back(p_to_e);
    }

    return crt(x, h_i);
}

// compute discrete log for a group of prime power order
// https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm
mpz_class discrete_log_prime_power_parallel(mpz_class g, mpz_class b, mpz_class mod, mpz_class p, mpz_class e, int num_workers_bsgs)
{
    if (e == 1) {
        return BabyStepGiantStep::discrete_log_parallel(g, b, mod, p, num_workers_bsgs);
    }
    mpz_class result = 0;

    const mpz_class gamma = powerMod(g, powerMod(p, e-1, mod), mod);

    for (mpz_class k = 0; k < e; ++k) {
        const mpz_class exponent = powerMod(p, e-1-k, mod);
        const mpz_class h_k = powerMod((b * powerMod(modInversePrime(g, mod), result, mod)) % mod, exponent, mod);
        const mpz_class dk = BabyStepGiantStep::discrete_log_parallel(gamma, h_k, mod, p, num_workers_bsgs);

        result += (dk * powerMod(p, k, mod)) % mod;
    }
    return result;
}

mpz_class PohligHellman::discrete_log_parallel(const mpz_class g, const mpz_class b, const mpz_class p, const std::vector<std::pair<mpz_class, mpz_class>> &factors) {
    const mpz_class order = p - 1;

    std::vector<mpz_class> x;
    x.reserve(factors.size());

    std::vector<mpz_class> h_i;
    h_i.reserve(factors.size());

    auto worker_body = [&](auto factor) {
        mpz_class prime = factor.first;
        mpz_class e = factor.second;

        mpz_class p_to_e = powerMod(prime, e, p);
        mpz_class subgroup_order = order / p_to_e;

        mpz_class subgroup_gen = powerMod(g, subgroup_order, p);
        mpz_class h = powerMod(b, subgroup_order, p);

        mpz_class x_i = discrete_log_prime_power_parallel(subgroup_gen, h, p, prime, e, 4); // TODO num_workers_bsgs

        x.push_back(x_i);
        h_i.push_back(p_to_e);
    };

    std::vector<std::thread> threads;
    for (auto factor : factors) {
        threads.push_back(std::thread(worker_body, factor));
    }

    for (auto &thread : threads) {
        thread.join();
    }

    return crt(x, h_i);
}