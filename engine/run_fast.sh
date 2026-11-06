#!/bin/bash
# run_fast.sh ID [SHARD NSHARDS]: exhaustive search of type ID with both reductions in the fast setting
# (output orbits on the first 10 pairs, input checkpoints at levels 3-6 with 8 maps), politely; output in out2/
id=$1
if [ -n "$2" ]; then tag=t_$id.$2; sh="-shard $2 $3 -sdepth 8"; else tag=t_$id; sh=""; fi
[ -f out2/$tag.done ] && exit 0
s0=$(date +%s)
./polite.sh ./sym10cp -dfs -dedupe -symdepth 10 -cp 3 6 -nsin 8 $sh < W3/w_$id.txt > out2/$tag.out 2> out2/$tag.err
echo "$(grep '^done' out2/$tag.err) secs $(( $(date +%s)-s0 ))" > out2/$tag.done
