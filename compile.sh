#!/bin/bash
gcc -W -Wall -Wunused -Wuninitialized -mpopcnt -O9 *.c -o ./slaMEM -lm
cp ./slaMEM /scratch/fjdf/mems/slaMEM/
gcc -W -Wall -Wunused -Wuninitialized -mpopcnt -O9 *.c -o ./slaMEM_t -lm -DBENCHMARK=1
mv ./slaMEM_t /scratch/fjdf/mems/slaMEM/
gcc -W -Wall -Wunused -mpopcnt -g -ggdb *.c -o ./debug/slaMEM-debug -lm
