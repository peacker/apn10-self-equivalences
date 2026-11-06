#!/bin/bash
# cost of one symmetry type: cost_one.sh LINE (a line of types10.tsv); appends to costs.tsv
IFS=$'\t' read -r id p kind fa fb A B <<< "$1"
w=W/w_$id.txt
info=$(python3 ../types/sym_setup.py 10 --amat "$A" --bmat "$B" -o $w)
d=$(echo "$info" | sed 's/.*d = \([0-9]*\);.*/\1/'); reps=$(echo "$info" | sed 's/.*C(B): \([0-9]*\);.*/\1/')
piv=$(echo "$info" | sed 's/.*per level: //' | tr -d ' ')
k=$(timeout 600 ./sym10b -knuth 3000 -seed $id < $w 2>&1 | grep "^total")
nodes=$(echo "$k" | awk '{print $3}' | tr -d ,); lv=$(echo "$k" | awk '{print $NF}')
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$p" "$kind" "$fa" "$fb" "$d" "$reps" "$piv" "${nodes:-timeout}" "${lv:--}" >> costs.tsv
