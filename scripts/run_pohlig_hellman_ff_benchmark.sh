#!/bin/bash

# test configuration 1 1 1
for nw_poh in 1 2 4 8 16 32
do
    for f in input/pohlig_hellman_*; do
        echo "# Running with ${nw_poh} workers in the farm, bsgs pipeline workers: 1 1 1 input file: $f"
        LD_LIBRARY_PATH=./lib/ ./bin/pohlig_hellman_ff $f ${nw_poh} 1 1 1
    done
done

# test configuration 2 2 2
for nw_poh in 1 2 4 8 16 32
do
    for f in input/pohlig_hellman_*; do
        echo "# Running with ${nw_poh} workers in the farm, bsgs pipeline workers: 2 2 2 input file: $f"
        LD_LIBRARY_PATH=./lib/ ./bin/pohlig_hellman_ff $f ${nw_poh} 2 2 2
    done
done

# test configuration 4 4 4
for nw_poh in 1 2 4 8 16 32
do
    for f in input/pohlig_hellman_*; do
        echo "# Running with ${nw_poh} workers in the farm, bsgs pipeline workers: 4 4 4 input file: $f"
        LD_LIBRARY_PATH=./lib/ ./bin/pohlig_hellman_ff $f ${nw_poh} 4 4 4
    done
done

# test configuration 1 2 8
for nw_poh in 1 2 4 8 16 32
do
    for f in input/pohlig_hellman_*; do
        echo "# Running with ${nw_poh} workers in the farm, bsgs pipeline workers: 1 2 8 input file: $f"
        LD_LIBRARY_PATH=./lib/ ./bin/pohlig_hellman_ff $f ${nw_poh} 1 2 8
    done
done