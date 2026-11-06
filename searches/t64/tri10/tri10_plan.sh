#!/bin/bash
# Structured n = 10 search, one configuration after the other, 140 shards each (nice 19).
cd "$(dirname "$0")"   # build first: cc -O3 -o tri10 tri10.c
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run() { # base terms field
  for s in $(seq 0 139); do echo "$1 $2 $3 $s 140"; done | xargs -P 140 -L 1 nice -n 19 python3 tri10_job.py
  echo "finished base=$1 terms=$2 field=$3 $(date)" >> results/progress.txt
}
mkdir -p results
for t in 5 6 7 8 9 10; do run 3 $t 4; done
for t in 1 2 3 4; do run 3 $t 32; done
for t in 5 6 7 8 9; do run 9 $t 4; done
for t in 5 6 7 8 9; do run 33 $t 4; done
echo ALL DONE >> results/progress.txt
