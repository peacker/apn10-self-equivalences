#!/bin/bash
# run41t.sh R J K: type 41, shard R of 10 of the first run (depth 8) split again: sub-shard J of K at depth 16 (sym10cp2 =
# sym10.c with -pshard/-pdepth; same options as run_cfg.sh cp). The K sub-shards of R partition shard R exactly.
r=$1; j=$2; k=$3; tag=r$r.$j.$k; [ -f out2/$tag.done ] && exit 0
s0=$(date +%s)
LIMIT=250 ./polite.sh ./sym10cp2 -dfs -dedupe -symdepth 10 -cp 3 6 -nsin 8 -nsamp 4 -pshard $r 10 -pdepth 8 -shard $j $k -sdepth 16 < W3/w_41.txt > out2/$tag.out 2> out2/$tag.err
echo "rc $? $(grep "^done" out2/$tag.err) secs $(( $(date +%s)-s0 ))" > out2/$tag.done
