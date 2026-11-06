#!/bin/bash
# ea_job.sh TAG CHUNK NCHUNKS: EA check of chunk CHUNK of out/TAG.out.gz -> ea/TAG.CHUNK.txt (+ .done)
tag=$1; c=$2; k=$3; cd "$(dirname "$0")"; f=ea/$tag.$c
[ -f $f.done ] && exit 0
s0=$(date +%s)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 LIMIT=250 ./polite.sh python3 ea_leaves.py out/$tag.out.gz $c $k > $f.txt 2> $f.err
echo "rc $? secs $(( $(date +%s)-s0 ))" > $f.done
