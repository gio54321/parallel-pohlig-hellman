#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cstdint>
#include <gmpxx.h>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>

#include "pohlig_hellman.h"
#include "bsgs.h"
#include "discrete_utils.h"


/*
 * Sequential version of the Pohlig-Hellman algorithm on groups of prime power order
 */
mpz_class discrete_log_prime_power(mpz_class g, mpz_class b, mpz_class mod, mpz_class p, mpz_class e)
{
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

/*
 * Sequential version of the general Pohlig-Hellman algorithm
 */
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

/*
 * Parallel version of the Pohlig-Hellman algorithm on groups of prime power order
 * it just calls the parallel version of the baby-step giant-step algorithm at each iteration
 */
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

/*
 * Parallel version of the general Pohlig-Hellman algorithm
 * 
 * num_workers is the number of workers that will compute Pohlig-Hellman iterations
 * num_workers_bsgs is the number of workers for each baby-step giant-step computation
 * total parallelism degree is num_workers * num_workers_bsgs
 */
mpz_class PohligHellman::discrete_log_parallel(const mpz_class g, const mpz_class b, const mpz_class p, const std::vector<std::pair<mpz_class, mpz_class>> &factors, int num_workers, int num_workers_bsgs) {
    const mpz_class order = p - 1;

    std::vector<mpz_class> x(factors.size());
    std::vector<mpz_class> h_i(factors.size());

    // job queue and lock, each job is a tuple of (job_index, prime, exponent)
    std::queue<std::tuple<int, mpz_class, mpz_class>> jobs;
    std::mutex jobs_mutex;

    auto worker_body = [&]() {
        while (true) {
            mpz_class prime;
            mpz_class e;
            int job_index;

            {
                // take a job from the shared queue
                std::lock_guard<std::mutex> lock(jobs_mutex);
                if (jobs.empty()) {
                    return;
                }
                std::tie(job_index, prime, e) = jobs.front();
                jobs.pop();
            }

            // compute the Pohlig-Hellman iteration
            mpz_class p_to_e = powerMod(prime, e, p);
            mpz_class subgroup_order = order / p_to_e;

            mpz_class subgroup_gen = powerMod(g, subgroup_order, p);
            mpz_class h = powerMod(b, subgroup_order, p);

            mpz_class x_i = discrete_log_prime_power_parallel(subgroup_gen, h, p, prime, e, num_workers_bsgs);

            // store the result
            x[job_index] = x_i;
            h_i[job_index] = p_to_e;
        }
    };

    for (size_t i = 0; i < factors.size(); ++i) {
        jobs.push(std::make_tuple(i, factors[i].first, factors[i].second));
    }

    std::vector<std::thread> threads;
    for (int i=0; i < num_workers; ++i) {
        threads.push_back(std::thread(worker_body));
    }

    for (auto &thread : threads) {
        thread.join();
    }

    // finally combine the results using the Chinese Remainder Theorem
    return crt(x, h_i);
}