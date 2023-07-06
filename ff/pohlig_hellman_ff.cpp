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

#include "discrete_utils.h"
#include "oatable.hpp"

#include "dlog_tasks.hpp"

using namespace ff;


struct TableInitializationStage : ff_Map<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        parallel_for(0,t->table.size(),[&](const long i) { 
            t->table[i] = std::numeric_limits<uint64_t>::max();
		}, 1); // TODO nw
        return t;
    }
};

// TODO mention in the report that an interface like parallel_for with (start, end) could be useful
// for tasks like these
struct TableConstructionStage : ff_Map<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        const int num_workers = 1;
        const int task_granularity = 4;
        size_t chunk_size = t->m.get_ui() / num_workers / task_granularity;
        if (chunk_size == 0) {
            chunk_size = 1;
        }

        parallel_for(0,t->m.get_ui(), chunk_size, 1, [&](long start) {
            const long end = std::min(start + chunk_size, t->m.get_ui());
            // std::cout << "[construction] start: " << start << ", stop: " << end << std::endl;

            if (start == 0) {
                table_insert(t->table, 1, 0);
                // std::cout << "inserted 0" << std::endl;
                start += 1;
            }
            mpz_class val = powerMod(t->g, mpz_class(start), t->p);
            for (size_t i = start; i < end; ++i){
                table_insert(t->table, val, i);
                // std::cout << "inserted " << val << "i: " << i << std::endl;
                val = (val * t->g) % t->p;
            }
        }, num_workers); // TODO nw
        return t;
    }
};

struct TableSearchStage : ff_Map<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        const int num_workers = 1;
        const int task_granularity = 4;
        size_t chunk_size = t->m.get_ui() / num_workers / task_granularity;
        if (chunk_size == 0) {
            chunk_size = 1;
        }

        std::mutex result_lock;
        std::atomic_bool found = false;


        parallel_for(0,t->m.get_ui(), chunk_size, 1, [&](long start) {
            const long end = std::min(start + chunk_size, t->m.get_ui());
            // std::cout << "start: " << start << ", stop: " << end << std::endl;

            mpz_class gm = (powerMod(modInversePrime(t->g, t->p), t->m, t->p)) % t->p;
            mpz_class val = (t->b * powerMod(gm, start, t->p)) % t->p;

            for (size_t i = start; i < end; ++i) {
                if (found) {
                    return;
                }
                std::optional<mpz_class> result = table_search(t->table, val, t->g, t->p);
                if (result.has_value()) {
                    std::lock_guard<std::mutex> lock(result_lock);
                    if (found) {
                        return;
                    }
                    t->result = (i * t->m + result.value()) % t->order;
                    found = true;
                    return;
                }
                val = (val * gm) % t->p;
            }

        }, num_workers); // TODO nw
        return t;
    }
};

poh_task_t *parse_task() {
    poh_task_t *task = new poh_task_t;

    std::ifstream myfile;
    myfile.open("input/pohlig_hellman_inputs.txt");

    std::string line;

    std::getline(myfile, line);
    std::cout << line << std::endl;
    std::getline(myfile, line);
    std::string sg, sb, sp;
    std::istringstream(line) >> sg >> sb >> sp;
    task->g = mpz_class(sg);
    task->b = mpz_class(sb);
    task->p = mpz_class(sp);
    task->order = task->p - 1;

    // std::cout << "g: " << task->g << ", b: " << task->b << ", p: " << task->p << std::endl;

    std::getline(myfile, line);
    auto line_stream = std::istringstream(line);

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

    float load_factor = 4.0; // TODO
    bsgs_task->table = std::vector<std::atomic_uint64_t>((size_t)(bsgs_task->m.get_ui() * load_factor));

    return bsgs_task;
}

struct Emitter : ff_monode_t<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *input_task) {

        // forward tasks from the feedback channel
        if (input_task != nullptr) {
            // std::cout << "received task: " << input_task->g << std::endl;
            return input_task;
        }

        auto task = parse_task();

        // create instances of poh_pp_task_t and bsgs_task_t
        int i=0;
        for (auto factor : task->order_factorization) {
            poh_pp_task_t *t = new poh_pp_task_t;
            t->factor = factor.first;
            t->exponent = factor.second;
            t->iteration = 0;
            t->parent_task = task;
            t->result = mpz_class(0);
            t->job_index = i;


            mpz_class p_to_e = powerMod(t->factor, t->exponent, task->p);
            mpz_class subgroup_order = task->order / p_to_e;

            mpz_class subgroup_gen = powerMod(task->g, subgroup_order, task->p);
            mpz_class h = powerMod(task->b, subgroup_order, task->p);

            t->g = subgroup_gen;
            t->b = h;

            t->gamma = powerMod(t->g, powerMod(t->factor, t->exponent-1, task->p), task->p);

            bsgs_task_t *t2 = create_bsgs_task_from_poh_pp(t);

            // std::cout << "g: " << t2->g << ", b: " << t2->b << ", p: " << t2->p << ", m: " << t2->m << std::endl;
            // std::cout << "m: " << t2->m << std::endl;
            ff_send_out(t2);

            i++;

        }

        return GO_ON;
    }

    void eosnotify(ssize_t t) {
        // when we receive EOS, send EOS to all workers of farm
        broadcast_task(EOS);
        return;
    }
};

struct Collector : ff_monode_t<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        // for (int i=0;i<t->table.size();++i) {
        //     std::cout << t->table[i] << " ";
        // }
        // std::cout << "received " << t->result;
        // std::cout << " g: " << t->g << ", b: " << t->b << ", p: " << t->p << ", m: " << t->m << std::endl;
        assert(powerMod(t->g, t->result, t->p) == t->b);

        if (t->parent_task != nullptr) {
            t->parent_task->result += t->result * powerMod(t->parent_task->factor, t->parent_task->iteration, t->p);
            t->parent_task->result %= t->p;

            t->parent_task->iteration += 1;
            if (t->parent_task->iteration == t->parent_task->exponent) {
                // std::cout << "result: " << t->parent_task->result << std::endl;
                t->parent_task->parent_task->x[t->parent_task->job_index] = t->parent_task->result;
                t->parent_task->parent_task->h_i[t->parent_task->job_index] = powerMod(t->parent_task->factor, t->parent_task->exponent, t->p);

                t->parent_task->parent_task->completed_subgroups += 1;
                // std::cout << "completed subgroups: " << t->parent_task->parent_task->completed_subgroups << std::endl;
                // std::cout << "order factorization size: " << t->parent_task->parent_task->order_factorization.size() << std::endl;
                if (t->parent_task->parent_task->completed_subgroups == t->parent_task->parent_task->order_factorization.size()) {
                    // std::cout << "size of x: " << t->parent_task->parent_task->x.size() << std::endl;
                    // for(int i=0;i<t->parent_task->parent_task->x.size();++i) {
                    //     std::cout << "x: " << t->parent_task->parent_task->x[i] << ", h_i: " << t->parent_task->parent_task->h_i[i] << std::endl;
                    // }
                    mpz_class result = crt(t->parent_task->parent_task->x, t->parent_task->parent_task->h_i);
                    std::cout << "general P-H result: " << result << std::endl;
                    assert(powerMod(t->parent_task->parent_task->g, result, t->parent_task->parent_task->p) == t->parent_task->parent_task->b);

                    delete t->parent_task->parent_task;
                    delete t->parent_task;
                    delete t;
                    // send EOS to the emitter, terminate
                    ff_send_out_to(EOS, 0);
                    return EOS;
                }

                delete t->parent_task;
            } else {
                bsgs_task_t *t2 = create_bsgs_task_from_poh_pp(t->parent_task);
                // std::cout << "sending back to bsgs ";
                // std::cout << "g: " << t2->g << ", b: " << t2->b << ", p: " << t2->p << ", m: " << t2->m << std::endl;
                ff_send_out_to(t2, 0);
            }
        }

        delete t;
        return GO_ON;
    }
};

int main(int argc, char * argv[]) {


    const int n_farms = 3; // TODO

    std::vector<ff_node *> bsgs_instances;
    for (int i=0; i<n_farms; ++i) {
        auto worker = new ff_Pipe<bsgs_task_t>(
        make_unique<TableInitializationStage>(), // TODO remember to add nw
        make_unique<TableConstructionStage>(),
        make_unique<TableSearchStage>());
        bsgs_instances.push_back(worker);
    }

    Emitter e;
    Collector c;

    ff_farm bsgs_farm;
    bsgs_farm.add_emitter(&e);
    bsgs_farm.add_workers(bsgs_instances);
    bsgs_farm.add_collector(&c);
    bsgs_farm.cleanup_workers();

    bsgs_farm.wrap_around();

    if (bsgs_farm.run_and_wait_end()<0) {
        error("running pipe");
        return -1;
    }

    return 0;
}