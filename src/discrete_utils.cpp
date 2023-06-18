#include "discrete_utils.h"

mpz_class powerMod(mpz_class a, mpz_class b, mpz_class m)
{
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

mpz_class gcdExtended(mpz_class a, mpz_class b, mpz_class* x, mpz_class* y)
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

// Function to find modulo inverse of a
mpz_class modInverse(mpz_class A, mpz_class M)
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

// Function to find inverse of a modulo prime p
mpz_class modInversePrime(mpz_class A, mpz_class p)
{
    return powerMod(A, p - 2, p);
}

/*
 * Chinese Remainder Theorem
 *
 * coeffs: vector of coefficients
 * moduli: vector of moduli, assumed to be coprime
 *
 * returns: x such that x = coeffs[i] (mod moduli[i]) for all i
 */
mpz_class crt(std::vector<mpz_class> &coeffs, std::vector<mpz_class> &moduli) {
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


uint64_t mulMod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t result = 0; 
    a = a % m;
    while (b > 0) {
        // If b is odd, add 'a' to result
        result = (result + a*(b%2)) % m;
 
        // Multiply 'a' with 2
        a = (a * 2) % m;
 
        // Divide b by 2
        b /= 2;
    }
 
    // Return result
    return result % m;
}


uint64_t powerMod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t res = 1;
    a = a % m;
    while(b>0)
    {
        if(b%2 == 1)
            res = mulMod(res, a, m);
        a = mulMod(a, a, m);
        b = b >> 1;
    }
    return res % m;
}

uint64_t modInversePrime(uint64_t A, uint64_t p) {
    return powerMod(A, p - 2, p);
}

extern "C" uint64_t mulModAsm(uint64_t a, uint64_t b, uint64_t m);
asm(R"(
.globl mulModAsm
mulModAsm:
.cfi_startproc
    movq %rdx, %r10 
    movq %rdi, %rax
    mulq %rsi
    divq %r10
    movq %rdx, %rax
    ret
.cfi_endproc
)");