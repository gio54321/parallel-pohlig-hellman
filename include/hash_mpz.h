#ifndef HASH_MPZ_H_
#define HASH_MPZ_H_

#include <cstddef>
#include <gmpxx.h>

template<> struct std::hash<mpz_srcptr> {
    size_t operator()(const mpz_srcptr x) const;
};

template<> struct std::hash<mpz_t> {
    size_t operator()(const mpz_t x) const;
};

template<> struct std::hash<mpz_class> {
    size_t operator()(const mpz_class &x) const;
};

#endif /* HASH_MPZ_H_ */