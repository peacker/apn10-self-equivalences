#!/bin/bash
# run24.sh GE R K: type 24 exhaustive search, Gold exponent GE (3 or 5), part R of K (K = 2048 in the paper run,
# R = 0..2047, so 4,096 jobs), politely; leaves in out/gGE.R.lut, summary in out/gGE.R.out, exit code and time in
# out/gGE.R.done. Reconstructed from the job outputs (the original script lived on the server; same options:
# env GE, PART = R/K, MINKEY, OUT). Build: cc -O3 -DM=5 -o full24g c/full24g.c
# Paper run: for ge in 3 5; do for r in $(seq 0 2047); do echo "$ge $r 2048"; done; done | xargs -P 198 -L 1 ./run24.sh
ge=$1; r=$2; k=$3; tag=g$ge.$r; mkdir -p out; [ -f out/$tag.done ] && exit 0
s0=$(date +%s)
GE=$ge PART=$r/$k MINKEY=1 OUT=out/$tag.lut LIMIT=250 ../../engine/polite.sh ./full24g > out/$tag.out 2> out/$tag.err
echo "rc $? secs $(( $(date +%s)-s0 ))" > out/$tag.done
