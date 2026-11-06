#!/bin/bash
# run73a.sh J: type 73 rerun printing every leaf (sym10f73 -all; same search as the 2026-11-02 run): shard J of 3072
# at depth 10. Shards J = r + 192 i (i = 0..15) partition piece r of the old 192-piece run exactly.
j=$1; tag=a$j; [ -f out/$tag.done ] && exit 0
[ -f claimed/$tag ] && exit 0
s0=$(date +%s)
LIMIT=250 ./polite.sh ./sym10f73 -dfs -dedupe -all -symdepth 45 -cp 3 10 -nsin 24 -nsamp 12 -shard $j 3072 -sdepth 10 < w_73.txt > out/$tag.out 2> out/$tag.err
echo "rc $? $(grep "^done" out/$tag.err) secs $(( $(date +%s)-s0 ))" > out/$tag.done
