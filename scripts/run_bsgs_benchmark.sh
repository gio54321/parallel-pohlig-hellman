#!/bin/bash

for i in {1..64}
do
    echo "# Running with $i workers"
    for f in input/bsgs_*; do
        echo "# Running with $f"
        LD_LIBRARY_PATH=./lib/ ./bin/test_bsgs $f $i 3
    done
done