#!/bin/bash
# exhaustive search of one symmetry type (or one shard of it): run_type.sh ID [SHARD NSHARDS]
id=$1; [ -f STOP_SYM10 ] && exit 0
if [ -n "$2" ]; then tag=t_$id.$2; extra="-shard $2 $3 -sdepth 10"; else tag=t_$id; extra=""; fi
[ -f out/$tag.done ] && exit 0
s0=$(date +%s)
./sym10b -dfs -dedupe $extra < W/w_$id.txt > out/$tag.out 2> out/$tag.err
echo "$(tail -1 out/$tag.err) secs $(( $(date +%s)-s0 ))" > out/$tag.done
