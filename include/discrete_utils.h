#ifndef DISCRETE_UTILS_H
#define DISCRETE_UTILS_H

#include <gmpxx.h>
#include <cstdint>
#include <vector>
#include <gmpxx.h>

/*
 * Modular binary exponentiation, it computes a^b mod m
 * in O(log(b)) time
 */
inline mpz_class powerMod(mpz_class a, mpz_class b, mpz_class m)
{
    if (b == 0) {
        return mpz_class(1);
    }
    mpz_class res = 1;
    a = a % m;
    while(b>0)
    {
        if(b%2 == 1)
            res = (res * a) % m;
        a = (a * a) % m;
        b = b >> 1;
    }
    return res % m;
}

/*
 * Compute gcd(a, b) and x, y such that ax + by = gcd(a, b)
 */
inline mpz_class gcdExtended(mpz_class a, mpz_class b, mpz_class* x, mpz_class* y)
{
    // Base Case
    if (a == 0) {
        *x = 0, *y = 1;
        return b;
    }
    // To store results of recursive call
    mpz_class x1, y1;
    mpz_class gcd = gcdExtended(b % a, a, &x1, &y1);
    // Update x and y using results of recursive
    // call
    *x = y1 - (b / a) * x1;
    *y = x1;
    return gcd;
}

/*
 * Compute inverse of a modulo m
 */
inline mpz_class modInverse(mpz_class A, mpz_class M)
{
    mpz_class x, y;
    mpz_class g = gcdExtended(A, M, &x, &y);
    if (g != 1)
        return 0;
    else {
        mpz_class res = (x % M + M) % M;
        return res;
    }
}

/*
 * Compute inverse of a modulo p, where p is prime
 * this case is faster than the general case, because we can use Fermat's little theorem
 * to compute a^(p-2) mod p, which is the inverse of a mod p
 */
inline mpz_class modInversePrime(mpz_class A, mpz_class p)
{
    return powerMod(A, p - 2, p);
}

/*
 * compute the solutions following the chinese Remainder Theorem
 *
 * coeffs: vector of coefficients
 * moduli: vector of moduli, assumed to be coprime
 *
 * returns: x such that x = coeffs[i] (mod moduli[i]) for all i
 */
inline mpz_class crt(std::vector<mpz_class> &coeffs, std::vector<mpz_class> &moduli) {
    mpz_class prod = 1;
    for (size_t i = 0; i < moduli.size(); i++) {
        prod *= moduli[i];
    }
    mpz_class result = 0;
    for (size_t i = 0; i < moduli.size(); i++) {
        mpz_class pp = prod / moduli[i];
        mpz_class inv = modInverse(pp, moduli[i]);
        result = (result + coeffs[i] * pp * inv) % prod;
    }
    return result % prod;
}

/*
 * Define the hash() function for mpz_class
 * since hashing is used only after doing modular exponentiation, we can just use
 * the least significant 64 bits as the hash value
 */
template<> struct std::hash<mpz_class> {
    size_t operator()(const mpz_class &x) const;
};

inline size_t std::hash<mpz_class>::operator()(const mpz_class &x) const {
    return (size_t) x.get_ui();
}


#endif