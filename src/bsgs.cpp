#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cstdint>
#include <gmpxx.h>

#include "bsgs.h"
#include "hash_mpz.h"
#include "discrete_utils.h"


BabyStepGiantStep::BabyStepGiantStep()
{
}

BabyStepGiantStep::~BabyStepGiantStep()
{
}


mpz_class BabyStepGiantStep::discrete_log(mpz_class g, mpz_class b, mpz_class p)
{
    mpz_class m = sqrt(p) + 1;
    std::unordered_map<mpz_class, mpz_class> table = {}; 

    for (mpz_class i = 0; i < m; ++i) {
        mpz_class val = powerMod(g, i, p);
        table.insert({val, i});
    }

    mpz_class gm = (powerMod(modInversePrime(g, p), m, p)) % p;
    mpz_class val = b;

    for (mpz_class i = 0; i < m; ++i) {
        if (table.find(val) != table.end()) {
            return (i * m + table[val]) % p;
        }
        val = (val * gm) % p;
    }
    return 0;
}

mpz_class BabyStepGiantStep::discrete_log_2(mpz_class g, mpz_class b, mpz_class p)
{
    mpz_class m = sqrt(p) + 1;
    std::unordered_map<mpz_class, mpz_class> table = {}; 
    table.reserve(m.get_ui());

    mpz_class val = g;
    table.insert({1, 0});
    for (mpz_class i = 1; i < m; ++i) {
        table.insert({val, i});
        val = (val * g) % p;
    }

    mpz_class gm = (powerMod(modInversePrime(g, p), m, p)) % p;
    val = b;

    for (mpz_class i = 0; i < m; ++i) {
        if (table.find(val) != table.end()) {
            return (i * m + table[val]) % p;
        }
        val = (val * gm) % p;
    }
    return mpz_class(0);
}

uint64_t BabyStepGiantStep::discrete_log_native(uint64_t g, uint64_t b, uint64_t p)
{
    mpz_class mz = (sqrt(mpz_class(p)) + 1);
    uint64_t m = mz.get_ui();

    // TODO check if this has the right bounds
    std::unordered_map<uint64_t, uint32_t> table = {}; 
    table.reserve(m);

    uint64_t val = g;
    table.insert({1, 0});
    for (uint32_t i = 1; i < m; ++i) {
        table.insert({val, i});
        val = mulMod(val, g, p);
    }

    uint64_t gm = (powerMod(modInversePrime(g, p), m, p)) % p;
    val = b;

    for (uint32_t i = 0; i < m; ++i) {
        if (table.find(val) != table.end()) {
            return (mulMod(i, m, p) + table[val]) % p;
        }
        val = mulMod(val, gm, p);
    }
    return 0; 
}

uint64_t BabyStepGiantStep::discrete_log_native_asm(uint64_t g, uint64_t b, uint64_t p)
{
    mpz_class mz = (sqrt(mpz_class(p)) + 1);
    uint64_t m = mz.get_ui();

    // TODO check if this has the right bounds
    std::unordered_map<uint64_t, uint32_t> table = {}; 
    table.reserve(m);

    uint64_t val = g;
    table.insert({1, 0});
    for (uint32_t i = 1; i < m; ++i) {
        table.insert({val, i});
        val = mulModAsm(val, g, p);
    }

    uint64_t gm = (powerMod(modInversePrime(g, p), m, p)) % p;
    val = b;

    for (uint32_t i = 0; i < m; ++i) {
        if (table.find(val) != table.end()) {
            return (mulModAsm(i, m, p) + table[val]) % p;
        }
        val = mulModAsm(val, gm, p);
    }
    return 0; 
}