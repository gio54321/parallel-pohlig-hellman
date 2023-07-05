#include <ff/ff.hpp>
#include <ff/map.hpp>
#include <iostream>
#include <atomic>
#include <vector>
#include <optional>
#include <algorithm>


#include <gmpxx.h>

#include "discrete_utils.h"
#include "oatable.hpp"

using namespace ff;

struct bsgs_task_t {
    int i;
    mpz_class g;
    mpz_class b;
    mpz_class p;
    mpz_class order;

    mpz_class m;
    mpz_class result;

    std::vector<std::atomic_uint64_t> table;
};


struct tableInitializationStage : ff_Map<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        parallel_for(0,t->table.size(),[&](const long i) { 
            t->table[i] = std::numeric_limits<uint64_t>::max();
		}, 2); // TODO nw
        return t;
    }
};

// TODO mention in the report that an interface like parallel_for with (start, end) could be useful
// for tasks like these
struct tableConstructionStage : ff_Map<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        const int num_workers = 2;
        const int task_granularity = 4;
        size_t chunk_size = t->m.get_ui() / num_workers / task_granularity;

        parallel_for(0,t->m.get_ui(), chunk_size, 1, [&](long start) {
            const long end = std::min(start + chunk_size, t->m.get_ui()-1);
            //std::cout << "start: " << start << ", stop: " << end << std::endl;

            if (start == 0) {
                table_insert(t->table, 0, 1);
                start += 1;
            }
            mpz_class val = powerMod(t->g, mpz_class(start), t->p);
            for (size_t i = start; i < end; ++i){
                table_insert(t->table, val, i);
                val = (val * t->g) % t->p;
            }
        }, num_workers); // TODO nw
        return t;
    }
};

struct tableSearchStage : ff_Map<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        const int num_workers = 2;
        const int task_granularity = 4;
        size_t chunk_size = t->m.get_ui() / num_workers / task_granularity;

        parallel_for(0,t->m.get_ui(), chunk_size, 1, [&](long start) {
            const long end = std::min(start + chunk_size, t->table.size());
            // std::cout << "start: " << start << ", stop: " << end << std::endl;

            mpz_class gm = (powerMod(modInversePrime(t->g, t->p), t->m, t->p)) % t->p;
            mpz_class val = (t->b * powerMod(gm, start, t->p)) % t->p;

            for (size_t i = start; i < end; ++i) {
                std::optional<mpz_class> result = table_search(t->table, val, t->g, t->p);
                if (result.has_value()) {
                    t->result = (i * t->m + result.value()) % t->order;
                    return;
                }
                val = (val * gm) % t->p;
            }

        }, num_workers); // TODO nw
        return t;
    }
};

struct emitter : ff_node_t<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        for (int i=0;i<2;++i) {
            bsgs_task_t *t = new bsgs_task_t;
            // 1862809527 2976686452 4225447451
            t->g = mpz_class(1862809527);
            t->b = mpz_class(2976686452);
            t->p = mpz_class(4225447451);
            t->m = sqrt(t->p-1) + 1;
            t->order = t->p - 1;
            t->result = mpz_class(-1);
            
            float load_factor = 4.0; // TODO
            t->table = std::vector<std::atomic_uint64_t>((size_t)(t->m.get_ui() * load_factor));
            t->i = i;
            ff_send_out(t);
        }
        return EOS;
    }
};

struct collector : ff_node_t<bsgs_task_t> {
    bsgs_task_t *svc(bsgs_task_t *t) {
        // for (int i=0;i<t->table.size();++i) {
        //     std::cout << t->table[i] << " ";
        // }
        std::cout << "received " << t->result << std::endl;
        assert(powerMod(t->g, t->result, t->p) == t->b);
        delete t;
        return GO_ON;
    }
};


int main(int argc, char * argv[]) {
    int nworkers = 3;
    long number = 111;

    emitter emitter;
    tableInitializationStage tableInitializationStage;
    tableConstructionStage tableConstructionStage;
    tableSearchStage tableSearchStage;
    collector collector;

    ff_Pipe<> pipe(
        emitter,
        tableInitializationStage,
        tableConstructionStage,
        tableSearchStage,
        collector
    );

    if (pipe.run_and_wait_end()<0) {
        error("running pipe");
        return -1;
    }

    return 0;
}