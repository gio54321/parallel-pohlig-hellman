#!/bin/bash

for nw in 1 2 4 8 16 32
do
    for i in 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10
    do
        for f in input/bsgs_48*; do
            echo "# Running with ${nw} workers: $f, load factor $i"
            LD_LIBRARY_PATH=./lib/ ./bin/test_bsgs $f ${nw} $i 
        done
    done
done