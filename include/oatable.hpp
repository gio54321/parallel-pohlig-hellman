#ifndef OATABLE_H
#define OATABLE_H

/*
 * Hash table implementation using addressing with quadratic probing
 *
 * The table is a vector of atomic uint64_t or uint64_t values.
 * The table has to be initialized with all values set to std::numeric_limits<uint64_t>::max().
 */


#include <atomic>
#include <vector>
#include <optional>
#include <gmpxx.h>

#include "discrete_utils.h"

/*
 * Insert value into the table (atomic version)
 *
 * table: the table
 * pos: position of the value in the table
 * val: value to be inserted
 */
inline void table_insert(std::vector<std::atomic_uint64_t>& table, mpz_class pos, uint64_t val)
{
    // open addressing with quadratic probing
    const mpz_class hash = pos % table.size();
    const uint64_t table_index = hash.get_ui();

    bool exchanged = false;
    size_t offset = 0;
    while (!exchanged) {
        uint64_t expected = std::numeric_limits<uint64_t>::max();
        exchanged = table[(table_index + offset * offset) % table.size()].compare_exchange_strong(expected, val);
        offset++;
    }
}

/*
 * Search for value in the table (atomic version)
 *
 * table: the table
 * value: value to be searched
 * g: generator of the group
 * p: prime modulus of the group
 *
 * returns: the value of i if some g^i=value has been found in the table, std::nullopt otherwise
 */
inline std::optional<mpz_class> table_search(std::vector<std::atomic_uint64_t>& table, mpz_class value, mpz_class g, mpz_class p)
{
    // open addressing with quadratic probing
    const mpz_class hash = value % table.size();
    const uint64_t index = hash.get_ui();

    size_t offset = 0;
    uint64_t table_value = table[(index + offset * offset) % table.size()];

    while (table_value != std::numeric_limits<uint64_t>::max()) {
        mpz_class guess = powerMod(g, mpz_class(table_value), p);
        if (guess == value) {
            // found value in the table
            return table_value;

        }
        offset++;
        table_value = table[(index + offset * offset) % table.size()];
    }
    return std::nullopt;
}


/*
 * Insert value into the table (non-atomic version)
 *
 * table: the table
 * pos: position of the value in the table
 * val: value to be inserted
 */
inline void table_insert(std::vector<std::uint64_t>& table, mpz_class pos, uint64_t val)
{
    // open addressing with quadratic probing
    const mpz_class hash = pos % table.size();
    const uint64_t table_index = hash.get_ui();

    size_t offset = 0;
    for (;;) {
        if(table[(table_index + offset * offset) % table.size()] == std::numeric_limits<uint64_t>::max()) {
            table[(table_index + offset * offset) % table.size()] = val;
            break;
        }
        offset++;
    }
}

/*
 * Search for value in the table (non-atomic version)
 *
 * table: the table
 * value: value to be searched
 * g: generator of the group
 * p: prime modulus of the group
 *
 * returns: the value of i if some g^i=value has been found in the table, std::nullopt otherwise
 */
inline std::optional<mpz_class> table_search(std::vector<std::uint64_t>& table, mpz_class value, mpz_class g, mpz_class p)
{
    // open addressing with quadratic probing
    const mpz_class hash = value % table.size();
    const uint64_t index = hash.get_ui();

    size_t offset = 0;
    uint64_t table_value = table[(index + offset * offset) % table.size()];

    while (table_value != std::numeric_limits<uint64_t>::max()) {
        mpz_class guess = powerMod(g, mpz_class(table_value), p);
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