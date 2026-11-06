#!/bin/bash
# plain forced DFS for n = 8 type 29, seeds 5 and 6, 12 shards each (validation of -gbasis + -symdepth)
cd "$(dirname "$0")/.."
for s in 5 6; do for i in $(seq 0 11); do echo "$s $i"; done; done | xargs -P 12 -n 2 sh -c 'nice -n 19 ./sym10f -dfs -dedupe -shard $1 12 -sdepth 6 -force val6/f8_$0.txt < val6/w8_29.txt > val6/plain_$0.$1.txt 2> val6/plain_$0.$1.err'
