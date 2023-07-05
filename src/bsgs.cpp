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


mpz_class BabyStepGiantStep::discrete_log(mpz_class g, mpz_class b, mpz_class p)
{
    return BabyStepGiantStep::discrete_log(g, b, p, p - 1);
}

mpz_class BabyStepGiantStep::discrete_log(mpz_class g, mpz_class b, mpz_class p, mpz_class order)
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
            // std::cout << "i: " << i << ", m: " << m << ", val: " << val << ", table[val]: " << table[val] << std::endl;
            return (i * m + table[val]) % order;
        }
        val = (val * gm) % p;
    }
    return mpz_class(0);
}




mpz_class BabyStepGiantStep::discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers) {
    return BabyStepGiantStep::discrete_log_parallel(g, b, p, order, num_workers, BabyStepGiantStep::default_load_factor);
}

mpz_class BabyStepGiantStep::discrete_log_parallel(mpz_class g, mpz_class b, mpz_class p, mpz_class order, int num_workers, float load_factor)
{
    if (order < num_workers) {
        return BabyStepGiantStep::discrete_log(g, b, p, order);
    }

    const mpz_class m = sqrt(order-1) + 1;
    const mpz_class slice_size = m / num_workers;

    std::vector<std::atomic_uint64_t> table((size_t)(m.get_ui() * load_factor));

    size_t table_slice = table.size() / num_workers;

    std::vector<std::thread> threads(num_workers);

    std::barrier finished_initializing_table(num_workers);
    std::barrier finished_computing_table(num_workers);
    bool collision_found = false;

    std::mutex result_lock;
    mpz_class result = mpz_class(-1);

    std::mutex debug_lock;

    auto worker_body = [&](const int num_worker) {
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

        for (size_t i = table_slice_start; i < table_slice_end; ++i) {
            table[i] = std::numeric_limits<uint64_t>::max();
        }

        finished_initializing_table.arrive_and_wait();

        if (num_worker == 0) {
            start = 1;
            table_insert(table, 1, 0);
        }

        // std::cout << "Worker " << num_worker << " computing from " << start << " to " << end << std::endl;

        mpz_class val = powerMod(g, start, p);
        for (mpz_class i = start; i < end; ++i) {
            // {
            //     std::lock_guard<std::mutex> lock(debug_lock);
            //     std::cout << "Worker " << num_worker << " computing " << i << " with value " << val << std::endl;
            // }
            table_insert(table, val, i.get_ui());
            val = (val * g) % p;
        }

        finished_computing_table.arrive_and_wait();

        // for (size_t i = 0; i < table.size(); ++i) {
        //     std::cout << "table[" << i << "] = " << table[i] << std::endl;
        // }

        const mpz_class gm = (powerMod(modInversePrime(g, p), m, p)) % p;

        val = b;
        if (num_worker == 0) {
            start = 0;
        } else {
            val *= powerMod(gm, start, p);
            val %= p;
        }


        for (mpz_class i = start; i < end; ++i) {
            if (collision_found) {
                return;
            }

            // std::cout << "Worker " << num_worker << " computing " << i << " with value " << val << std::endl;
            std::optional<mpz_class> table_value = table_search(table, val, g, p);
            if (table_value.has_value()) {
                std::lock_guard<std::mutex> lock(result_lock);
                if (collision_found) {
                    return;
                }
                // std::cout << "Collision found by worker " << num_worker << " with i = " << i << std::endl;
                // std::cout <<"table_value: " << table_value.value() << ", val: " << val << std::endl;

                result = (i * m + table_value.value()) % order;
                collision_found = true;
                return;
            }
            val = (val * gm) % p;
        }
    };

    for (int i = 0; i < num_workers; ++i) {
        threads[i] = std::thread(worker_body, i);
    }

    for (int i = 0; i < num_workers; ++i) {
        threads[i].join();
    }
    return result;
}