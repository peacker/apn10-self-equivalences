#!/bin/bash
# run_cfg.sh ID SHARD NSHARDS MODE: MODE nocp = output orbits on 10 pairs only; cp = plus input checkpoints at
# levels 3-6 (8 maps), 4 sampled output maps per orbit computation. Politely; output in out3/
id=$1; i=$2; k=$3; mode=$4
tag=t_$id.$i; [ -f out3/$tag.done ] && exit 0
if [ "$mode" = cp ]; then opt="-symdepth 10 -cp 3 6 -nsin 8 -nsamp 4"; else opt="-symdepth 10 -nsamp 4"; fi
s0=$(date +%s)
./polite.sh ./sym10cp -dfs -dedupe $opt -shard $i $k -sdepth 8 -progress out3/$tag.prog < W3/w_$id.txt > out3/$tag.out 2> out3/$tag.err
echo "$(grep '^done' out3/$tag.err) secs $(( $(date +%s)-s0 ))" > out3/$tag.done
