#include <gmpxx.h>
#include <vector>
#include <atomic>


struct poh_task_t {
    mpz_class g;
    mpz_class b;
    mpz_class p;
    mpz_class order;
    std::vector<std::pair<mpz_class, mpz_class>> order_factorization;

    std::vector<mpz_class> x;
    std::vector<mpz_class> h_i;

    int completed_subgroups;
};

struct poh_pp_task_t {
    mpz_class g;
    mpz_class b;
    mpz_class factor;
    mpz_class exponent;
    mpz_class iteration;
    mpz_class gamma;

    mpz_class result;
    int job_index;

    poh_task_t *parent_task;
};

struct bsgs_task_t {
    mpz_class g;
    mpz_class b;
    mpz_class p;
    mpz_class order;

    mpz_class m;
    mpz_class result;

    std::vector<std::atomic_uint64_t> table;

    // pointer to parent Poligh-Hellman task
    // TODO is this a unique pointer?
    poh_pp_task_t *parent_task;
};