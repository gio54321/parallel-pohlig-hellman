#ifndef DISCRETE_UTILS_H
#define DISCRETE_UTILS_H

#include <gmpxx.h>
#include <cstdint>
#include <vector>
#include <gmpxx.h>

mpz_class powerMod(mpz_class a, mpz_class b, mpz_class m);
mpz_class gcdExtended(mpz_class a, mpz_class b, mpz_class* x, mpz_class* y);
mpz_class modInverse(mpz_class A, mpz_class M);
mpz_class modInversePrime(mpz_class A, mpz_class p);
mpz_class crt(std::vector<mpz_class> &coeffs, std::vector<mpz_class> &moduli);

// implement hash() for mpz_class
// since hashing is used only after doing modular exponentiation, we can just use the least significant 64 bits
// as the hash value
template<> struct std::hash<mpz_class> {
    size_t operator()(const mpz_class &x) const;
};


#endif