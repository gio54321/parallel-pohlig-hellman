#!/bin/bash

for nw_bsgs in 1 2 4 8
do

    for nw_poh in {0..32}
    do
        for f in input/pohlig_hellman_*; do
            echo "# Running with ${nw_poh} Pohlig-Hellman workers, ${nw_bsgs} bsgs workers: $f"
            LD_LIBRARY_PATH=./lib/ ./bin/test_pohlig_hellman $f ${nw_poh} ${nw_bsgs}
        done
    done
done