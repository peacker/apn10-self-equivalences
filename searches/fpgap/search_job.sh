#!/bin/bash
# search_job.sh TYPE [SHARD]: rerun of the original sym10b run (same options) printing every leaf (-all); gzip output
t=$1; s=$2; cd "$(dirname "$0")"
if [ -n "$s" ]; then tag=t_$t.$s; extra="-shard $s 16 -sdepth 10"; else tag=t_$t; extra=""; fi
[ -f out/$tag.done ] && exit 0
s0=$(date +%s)
LIMIT=250 ./polite.sh ./sym10f73 -dfs -dedupe -all $extra < W/w_$t.txt 2> out/$tag.err | gzip -1 > out/$tag.out.gz
echo "rc ${PIPESTATUS[0]} $(grep '^done' out/$tag.err) secs $(( $(date +%s)-s0 ))" > out/$tag.done
