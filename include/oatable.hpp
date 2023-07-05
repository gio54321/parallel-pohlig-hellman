#ifndef OATABLE_H
#define OATABLE_H

#include <atomic>
#include <vector>
#include <optional>
#include <gmpxx.h>

#include "discrete_utils.h"

inline void table_insert(std::vector<std::atomic_uint64_t>& table, mpz_class pos, uint64_t val)
{
    // open addressing with linear probing
    const mpz_class hash = pos % table.size();
    const uint64_t table_index = hash.get_ui();

    bool exchanged = false;
    size_t offset = 0;
    while (!exchanged) {
        uint64_t expected = std::numeric_limits<uint64_t>::max();
        exchanged = table[(table_index + offset * offset) % table.size()].compare_exchange_strong(expected, val);
        //std::cout << "table_index: " << table_index << ", expected: " << expected << ", exchanged: " << exchanged << std::endl;
        offset++;
    }
}

inline std::optional<mpz_class> table_search(std::vector<std::atomic_uint64_t>& table, mpz_class value, mpz_class g, mpz_class p)
{
    // open addressing with linear probing
    const mpz_class hash = value % table.size();
    const uint64_t index = hash.get_ui();

    size_t offset = 0;
    uint64_t table_value = table[(index + offset * offset) % table.size()];

    while (table_value != std::numeric_limits<uint64_t>::max()) {
        mpz_class guess = powerMod(g, mpz_class(table_value), p);
        // std::cout << "table_value " << table_value<<  " guess: " << guess << std::endl;

        if (guess == value) {
            // found value in the table
            return table_value;

        }
        offset++;
        table_value = table[(index + offset * offset) % table.size()];
    }
    return std::nullopt;
}

#endif