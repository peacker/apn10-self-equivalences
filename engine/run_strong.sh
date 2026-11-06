#!/bin/bash
# run_strong.sh ID SHARD NSHARDS: full symmetry setting (output orbits on all pairs, checkpoints at levels 3-10 with
# 24 input maps, 12 output samples), politely; output and progress in out4/
id=$1; i=$2; k=$3
tag=t_$id.$i; [ -f out4/$tag.done ] && exit 0
s0=$(date +%s)
./polite.sh ./sym10cp -dfs -dedupe -symdepth 45 -cp 3 10 -nsin 24 -nsamp 12 -shard $i $k -sdepth 8 -progress out4/$tag.prog < W3/w_$id.txt > out4/$tag.out 2> out4/$tag.err
echo "$(grep '^done' out4/$tag.err) secs $(( $(date +%s)-s0 ))" > out4/$tag.done
