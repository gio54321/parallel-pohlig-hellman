#include <ff/ff.hpp>
#include <ff/map.hpp>
#include <ff/farm.hpp>
#include <iostream>
#include <atomic>
#include <vector>
#include <optional>
#include <algorithm>
#include <sstream>
#include <mutex>
#include <gmpxx.h>

#include "utimer.h"
#include "discrete_utils.h"
#include "oatable.hpp"
#include "dlog_tasks.hpp"

using namespace ff;

#define DEFAULT_LOAD_FACTOR 4.5

/*
 * Table initialization stage, sets all values in the table to the sentinel value (max uint64_t)
 */
struct TableInitializationStage : ff_Map<bsgs_task_t> {
    int num_workers;
    TableInitializationStage(int num_workers):num_workers(num_workers) {}

    bsgs_task_t *svc(bsgs_task_t *t) {
        parallel_for(0,t->table.size(),[&](const long i) { 
            t->table[i] = std::numeric_limits<uint64_t>::max();
		}, num_workers);
        return t;
    }
};

/*
 * Table construction stage, inserts g^i for i in [0, m) into the table
 */
struct TableConstructionStage : ff_Map<bsgs_task_t> {
    int num_workers;
    TableConstructionStage(int num_workers):num_workers(num_workers) {}

    bsgs_task_t *svc(bsgs_task_t *t) {
        const int task_granularity = 4;
        size_t chunk_size = t->m.get_ui() / num_workers / task_granularity;
        if (chunk_size == 0) {
            chunk_size = 1;
        }

        parallel_for(0,t->m.get_ui(), chunk_size, 1, [&](long start) {
            const long end = std::min(start + chunk_size, t->m.get_ui());

            if (start == 0) {
                table_insert(t->table, 1, 0);
                start += 1;
            }
            mpz_class val = powerMod(t->g, mpz_class(start), t->p);
            for (size_t i = start; i < end; ++i){
                table_insert(t->table, val, i);
                val = (val * t->g) % t->p;
            }
        }, num_workers);
        return t;
    }
};

/*
 * Table search stage, searches for b * g^(-im) for i in [0, m) in the table
 */
struct TableSearchStage : ff_Map<bsgs_task_t> {
    int num_workers;
    TableSearchStage(int num_workers):num_workers(num_workers) {}

    bsgs_task_t *svc(bsgs_task_t *t) {
        const int task_granularity = 4;
        size_t chunk_size = t->m.get_ui() / num_workers / task_granularity;
        if (chunk_size == 0) {
            chunk_size = 1;
        }

        std::mutex result_lock;
        std::atomic_bool found = false;

        parallel_for(0,t->m.get_ui(), chunk_size, 1, [&](long start) {
            const long end = std::min(start + chunk_size, t->m.get_ui());

            mpz_class gm = (powerMod(modInversePrime(t->g, t->p), t->m, t->p)) % t->p;
            mpz_class val = (t->b * powerMod(gm, start, t->p)) % t->p;

            for (size_t i = start; i < end; ++i) {
                if (found) {
                    return;
                }
                std::optional<mpz_class> result = table_search(t->table, val, t->g, t->p);
                if (result.has_value()) {
                    // result found, set the result and return
                    std::lock_guard<std::mutex> lock(result_lock);

                    // if the result has already been found, then return
                    if (found) {
                        return;
                    }

                    t->result = (i * t->m + result.value()) % t->order;
                    found = true;
                    return;
                }
                val = (val * gm) % t->p;
            }

        }, num_workers);
        return t;
    }
};


bsgs_task_t *create_bsgs_task_from_poh_pp(poh_pp_task_t *t) {
    bsgs_task_t *bsgs_task = new bsgs_task_t;

    const mpz_class exponent = powerMod(t->factor, t->exponent-1-t->iteration, t->parent_task->p);
    const mpz_class h_k = powerMod((t->b * powerMod(modInversePrime(t->g, t->parent_task->p), t->result, t->parent_task->p)) % t->parent_task->p, exponent, t->parent_task->p);

    bsgs_task->g = t->gamma;
    bsgs_task->b = h_k;
    bsgs_task->p = t->parent_task->p;
    bsgs_task->order = t->factor;
    bsgs_task->m = sqrt(bsgs_task->order-1) + 1;
    bsgs_task->result = mpz_class(-1);
    bsgs_task->parent_task = t;

    float load_factor = DEFAULT_LOAD_FACTOR;
    bsgs_task->table = std::vector<std::atomic_uint64_t>((size_t)(bsgs_task->m.get_ui() * load_factor));

    return bsgs_task;
}

struct Emitter : ff_monode_t<bsgs_task_t> {
    poh_task_t *pohlig_hellman_task;
    Emitter(poh_task_t *task):pohlig_hellman_task(task) {}

    bsgs_task_t *svc(bsgs_task_t *input_task) {

        // if there is a task from the feedback channel, just forward it to the farm
        if (input_task != nullptr) {
            return input_task;
        }

        // create all instances of poh_pp_task_t and bsgs_task_t
        int i=0;
        for (auto factor : pohlig_hellman_task->order_factorization) {
            poh_pp_task_t *t = new poh_pp_task_t;
            t->factor = factor.first;
            t->exponent = factor.second;
            t->iteration = 0;
            t->parent_task = pohlig_hellman_task;
            t->result = mpz_class(0);
            t->job_index = i;

            mpz_class p_to_e = powerMod(t->factor, t->exponent, pohlig_hellman_task->p);
            mpz_class subgroup_order = pohlig_hellman_task->order / p_to_e;
            mpz_class subgroup_gen = powerMod(pohlig_hellman_task->g, subgroup_order, pohlig_hellman_task->p);
            mpz_class h = powerMod(pohlig_hellman_task->b, subgroup_order, pohlig_hellman_task->p);
            t->g = subgroup_gen;
            t->b = h;
            t->gamma = powerMod(t->g, powerMod(t->factor, t->exponent-1, pohlig_hellman_task->p), pohlig_hellman_task->p);

            bsgs_task_t *t2 = create_bsgs_task_from_poh_pp(t);
            ff_send_out(t2);
            i++;

        }

        return GO_ON;
    }

    void eosnotify(ssize_t) {
        // when we receive EOS (from the feedback channel), send EOS to all workers of the farm
        broadcast_task(EOS);
        return;
    }
};

/*
 * Collector stage, collects the results of the bsgs tasks, computes the next iterations
 * and in the end recombines the results using crt
 */
struct Collector : ff_monode_t<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        // check that the bsgs discrete log has been computed correctly
        assert(powerMod(t->g, t->result, t->p) == t->b);


        if (t->parent_task != nullptr) {
            // update the parent task result
            t->parent_task->result += t->result * powerMod(t->parent_task->factor, t->parent_task->iteration, t->p);
            t->parent_task->result %= t->p;

            // increment the iteration number of parent task
            t->parent_task->iteration += 1;

            // check if this is the last iteration of the parent task
            if (t->parent_task->iteration == t->parent_task->exponent) {

                // put the resulting discrete log in the correct position in the general Pohlig
                // Hellman result vector
                t->parent_task->parent_task->x[t->parent_task->job_index] = t->parent_task->result;
                t->parent_task->parent_task->h_i[t->parent_task->job_index] = powerMod(t->parent_task->factor, t->parent_task->exponent, t->p);

                t->parent_task->parent_task->completed_subgroups += 1;

                // check if all subgroups have been computed
                if (t->parent_task->parent_task->completed_subgroups == t->parent_task->parent_task->order_factorization.size()) {
                    // combine the results using crt
                    mpz_class result = crt(t->parent_task->parent_task->x, t->parent_task->parent_task->h_i);
                    std::cout << "general P-H result: " << result << std::endl;

                    // check heck that we actually found the correct result for the whole task
                    assert(powerMod(t->parent_task->parent_task->g, result, t->parent_task->parent_task->p) == t->parent_task->parent_task->b);

                    // clean up previously allocated tasks
                    delete t->parent_task->parent_task;
                    delete t->parent_task;
                    delete t;

                    // send EOS to the emitter, terminate all workers
                    ff_send_out_to(EOS, 0);
                    return EOS;
                }

                delete t->parent_task;
            } else {
                // not the last iteration of the parent task, compute the next input parameters
                // for bsgs and send the task to the emitter
                bsgs_task_t *t2 = create_bsgs_task_from_poh_pp(t->parent_task);
                ff_send_out_to(t2, 0);
            }
        }

        delete t;
        return GO_ON;
    }
};

/*
 * Parse the input file and create the Pohlig-Hellman task
 */
poh_task_t *parse_task(char *filename) {
    poh_task_t *task = new poh_task_t;

    std::ifstream myfile;
    myfile.open(filename);

    if (!myfile.is_open()) {
        std::cout << "Error opening file" << std::endl;
    }

    std::string line;

    std::getline(myfile, line);
    std::string sg, sb, sp;
    std::istringstream(line) >> sg >> sb >> sp;
    task->g = mpz_class(sg);
    task->b = mpz_class(sb);
    task->p = mpz_class(sp);
    task->order = task->p - 1;

    std::getline(myfile, line);
    auto line_stream = std::istringstream(line);

    // read the factorization of the order of the group
    int size = 0;
    line_stream >> size;

    for (int i = 0; i < size; ++i) {
        std::string s_prime, s_e;
        line_stream >> s_prime >> s_e;
        mpz_class prime(s_prime);
        mpz_class e(s_e);
        task->order_factorization.push_back({prime, e});
    }

    task->x = std::vector<mpz_class>(task->order_factorization.size());
    task->h_i = std::vector<mpz_class>(task->order_factorization.size());
    task->completed_subgroups = 0;
    return task;
}

int main(int argc, char * argv[]) {

    if (argc != 6) {
        std::cout << "Usage: ./main <input_filename> <n_bsgs_farms> <nw_init> <nw_construction> <nw_search>" << std::endl;
        return 1;
    }

    char *filename = argv[1];
    int n_farms = atoi(argv[2]);
    int nw_init = atoi(argv[3]);
    int nw_construction = atoi(argv[4]);
    int nw_search = atoi(argv[5]);


    // create the bsgs workers
    std::vector<ff_node *> bsgs_instances;
    for (int i=0; i<n_farms; ++i) {
        auto worker = new ff_Pipe<bsgs_task_t>(
        make_unique<TableInitializationStage>(nw_init),
        make_unique<TableConstructionStage>(nw_construction),
        make_unique<TableSearchStage>(nw_search));
        bsgs_instances.push_back(worker);
    }

    Emitter e(parse_task(filename));
    Collector c;

    // create the farm
    ff_farm bsgs_farm;
    bsgs_farm.add_emitter(&e);
    bsgs_farm.add_workers(bsgs_instances);
    bsgs_farm.add_collector(&c);

    // set scheduling on demand, because the tasks are few
    // this could lead to better load balancing
    bsgs_farm.set_scheduling_ondemand();
    bsgs_farm.cleanup_workers();

    // create feedback channel
    bsgs_farm.wrap_around();

    long time_taken;
    {
        utimer timer("ff", &time_taken);

        if (bsgs_farm.run_and_wait_end()<0) {
            error("running pipe");
            return -1;
        }
    }

    std::cout << "Found result in " << time_taken << " microseconds" << std::endl;

    return 0;
}