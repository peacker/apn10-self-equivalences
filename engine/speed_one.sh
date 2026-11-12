#!/bin/bash
# speed_one.sh ID: nodes per second of the new engine on type ID (DFS limited to 3e8 nodes); appends to speeds.tsv
id=$1
s0=$(date +%s.%N)
r=$(./sym10c -dfs -limit 300000000 < W/w_$id.txt 2>&1 | grep '^done')
s1=$(date +%s.%N)
nodes=$(echo "$r" | sed 's/done: nodes \([0-9]*\),.*/\1/')
echo -e "$id\t$nodes\t$(echo "$s1 - $s0" | bc)" >> speeds.tsv
