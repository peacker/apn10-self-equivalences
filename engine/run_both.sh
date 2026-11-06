#!/bin/bash
# run_both.sh ID: exhaustive search of type ID with both symmetry reductions (politely), output in out/
id=$1; [ -f out/t_$id.done ] && exit 0
s0=$(date +%s)
./polite.sh ./sym10cp -dfs -dedupe -symdepth 45 -cp 3 10 < W3/w_$id.txt > out/t_$id.out 2> out/t_$id.err
echo "$(grep '^done' out/t_$id.err) secs $(( $(date +%s)-s0 ))" > out/t_$id.done
