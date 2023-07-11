#include <gmpxx.h>
#include <vector>
#include <atomic>

/*
 * Task strucure for the entire Pohlig-Hellman algorithm
 */
struct poh_task_t {
    // input parameters
    mpz_class g;
    mpz_class b;
    mpz_class p;
    mpz_class order;
    std::vector<std::pair<mpz_class, mpz_class>> order_factorization;

    // output vectors
    std::vector<mpz_class> x;
    std::vector<mpz_class> h_i;

    // number of completed subgroups, i.e., number of completed number of completed Pohlig-Hellman
    int completed_subgroups;
};

/*
 * Task structure for one computation of Pohlig-Hellman on a group of prime power order
*/
struct poh_pp_task_t {
    // input parameters, group order = factor^exponent
    mpz_class g;
    mpz_class b;
    mpz_class factor;
    mpz_class exponent;

    // computed gamma for the task
    mpz_class gamma;

    // iteration number
    mpz_class iteration;

    // result of the computation
    mpz_class result;

    // index of the job in the parent Pohlig-Hellman task
    int job_index;
    // pointer to parent Pohlig-Hellman task
    poh_task_t *parent_task;
};


/*
 * Task structure for one computation of Baby-Step Giant-Step algorithm
 */
struct bsgs_task_t {
    // input parameters
    mpz_class g;
    mpz_class b;
    mpz_class p;
    mpz_class order;

    // computed m for the task
    mpz_class m;

    // result of the computation
    mpz_class result;

    // open addressing hash table
    std::vector<std::atomic_uint64_t> table;

    // pointer to parent Pohlig-Hellman on prime power order groups task
    poh_pp_task_t *parent_task;
};