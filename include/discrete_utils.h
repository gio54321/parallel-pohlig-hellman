#include <gmpxx.h>
#include <cstdint>
#include <vector>
#include <gmpxx.h>

mpz_class powerMod(mpz_class a, mpz_class b, mpz_class m);
mpz_class gcdExtended(mpz_class a, mpz_class b, mpz_class* x, mpz_class* y);
mpz_class modInverse(mpz_class A, mpz_class M);
mpz_class modInversePrime(mpz_class A, mpz_class p);

mpz_class crt(std::vector<mpz_class> &coeffs, std::vector<mpz_class> &moduli);

// native implementations
uint64_t mulMod(uint64_t a, uint64_t b, uint64_t m);
extern "C" uint64_t mulModAsm(uint64_t a, uint64_t b, uint64_t m);

uint64_t powerMod(uint64_t a, uint64_t b, uint64_t m);
uint64_t modInversePrime(uint64_t A, uint64_t p);