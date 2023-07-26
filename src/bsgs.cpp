#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cstdint>
#include <gmpxx.h>
#include <thread>
#include <barrier>
#include <assert.h>
#include <atomic>
#include <map>
#include <mutex>

#include "bsgs.h"
#include "discrete_utils.h"
#include "oatable.hpp"

/*
 * Sequential version of the baby-step giant-step algorithm, using the open addressing hash table
 */
mpz_class BabyStepGiantStep::discrete_log(mpz_class g, mpz_class b, mpz_class p, mpz_class order)
{
    mpz_class m = sqrt(order-1) + 1;
    std::vector<std::uint64_t> table((size_t)(m.get_ui() * BabyStepGiantStep::default_load_factor));

    for (size_t i = 0; i < table.size(); ++i) {
        table[i] = std::numeric_limits<uint64_t>::max();
    }

    table_insert(table, 1, 0);
    mpz_class val = g;
    for (uint64_t i = 1; i < m; ++i) {
        table_insert(table, val, i);
        val = (val * g) % p;
    }

    mpz_class gm = (powerMod(modInversePrime(g, p), m, p)) % p;
    val = b;

    for (mpz_class i = 0; i < m; ++i) {
        std::optional<mpz_class> table_value = table_search(table, val, g, p);
        if (table_value.has_value()) {
            return (i * m + table_value.value()) % order;
        }
        val = (val * gm) % p;
    }
    return mpz_class(-1);
}


/*
 * Parallel version of the baby-step giant-step algorithm, using the open addressing hash table
 */
mpz_class BabyStepGiantStep::discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers, float load_factor)
{
    // if the order is less than the number of workers, there could be workers that don't do anything
    // so fall-back to the sequential version
    if (order < num_workers) {
        return BabyStepGiantStep::discrete_log(g, b, p, order);
    }

    const mpz_class m = sqrt(order-1) + 1;
    const mpz_class slice_size = m / num_workers;

    // initialize the tbale and compute the slice of the table that each worker will use
    std::vector<std::atomic_uint64_t> table((size_t)(m.get_ui() * load_factor));
    size_t table_slice = table.size() / num_workers;

    // initialize the vector of threads
    std::vector<std::thread> threads(num_workers);

    // initialize the barriers for the three stages of the algorithm, and the result flag and mutex
    std::barrier finished_initializing_table(num_workers);
    std::barrier finished_computing_table(num_workers);
    bool collision_found = false;
    std::mutex result_lock;

    mpz_class result = mpz_class(-1);

    auto worker_body = [&](const int num_worker) {
        // calculate start and end indices for the worker
        mpz_class start = num_worker * slice_size;
        mpz_class end = (num_worker + 1) * slice_size;

        if (num_worker == num_workers - 1) {
            end = m;
        }

        size_t table_slice_start = table_slice * num_worker;
        size_t table_slice_end = table_slice * (num_worker + 1);
        if (num_worker == num_workers - 1) {
            table_slice_end = table.size();
        }

        // initialize our slice of the table
        for (size_t i = table_slice_start; i < table_slice_end; ++i) {
            table[i] = std::numeric_limits<uint64_t>::max();
        }

        // block until all workers have finished initializing their slice of the table
        finished_initializing_table.arrive_and_wait();

        // insert the elements in the table
        if (num_worker == 0) {
            start = 1;
            table_insert(table, 1, 0);
        }

        mpz_class val = powerMod(g, start, p);
        for (mpz_class i = start; i < end; ++i) {
            table_insert(table, val, i.get_ui());
            val = (val * g) % p;
        }

        // block until all workers have finished inserting their elements in the table
        finished_computing_table.arrive_and_wait();

        // compute the search parameters
        const mpz_class gm = (powerMod(modInversePrime(g, p), m, p)) % p;

        val = b;
        if (num_worker == 0) {
            start = 0;
        } else {
            val *= powerMod(gm, start, p);
            val %= p;
        }


        for (mpz_class i = start; i < end; ++i) {
            // if a collision has been already found by another thread, stop the computation
            if (collision_found) {
                return;
            }

            // search for the value in the table
            std::optional<mpz_class> table_value = table_search(table, val, g, p);
            if (table_value.has_value()) {
                std::lock_guard<std::mutex> lock(result_lock);
                // it could be that two threads find a collision at the same time, so check again the flag
                if (collision_found) {
                    return;
                }

                // compute the result and set the collision flag
                result = (i * m + table_value.value()) % order;
                collision_found = true;
                return;
            }
            val = (val * gm) % p;
        }
    };

    // start the threads
    for (int i = 0; i < num_workers; ++i) {
        threads[i] = std::thread(worker_body, i);
    }

    // wait for the threads to finish
    for (int i = 0; i < num_workers; ++i) {
        threads[i].join();
    }
    return result;
}


/*
 * Old version of the baby-step giant-step algorithm, using std::unordered_map
 */
mpz_class BabyStepGiantStep::discrete_log_old(mpz_class g, mpz_class b, mpz_class p, mpz_class order)
{
    mpz_class m = sqrt(order-1) + 1;
    std::unordered_map<mpz_class, mpz_class> table = {};
    table.reserve(m.get_ui());
    table.insert({0, 1});

    mpz_class val = g;
    for (mpz_class i = 1; i < m; ++i) {
        table.insert({val, i});
        val = (val * g) % p;
    }

    mpz_class gm = (powerMod(modInversePrime(g, p), m, p)) % p;
    val = b;

    for (mpz_class i = 0; i < m; ++i) {
        if (table.find(val) != table.end()) {
            return (i * m + table[val]) % order;
        }
        val = (val * gm) % p;
    }
    return mpz_class(0);
}