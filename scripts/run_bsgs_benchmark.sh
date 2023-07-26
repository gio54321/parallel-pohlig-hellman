#!/bin/bash

for i in {0..64}
do
    for f in input/bsgs_*; do
        echo "# Running with $i workers: $f"
        LD_LIBRARY_PATH=./lib/ ./bin/test_bsgs $f $i 4.5
    done
done