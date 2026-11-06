#!/bin/bash
# knuth_sym.sh ID: regenerate W2/w_ID.txt (with the Z(B) basis) and estimate the tree with output symmetry
id=$1
line=$(grep -P "^$id\t" ../types/types10.tsv); A=$(echo "$line" | cut -f6); B=$(echo "$line" | cut -f7)
python3 ../types/sym_setup.py 10 --amat "$A" --bmat "$B" -o W2/w_$id.txt > /dev/null
k=$(timeout 1800 ./sym10s -knuth 2000 -seed $id -symdepth 45 < W2/w_$id.txt 2>&1 | grep "^total")
echo -e "$id\t$(echo "$k" | awk '{print $3}' | tr -d ,)\t$(echo "$k" | awk '{print $NF}')" >> knuth_sym.tsv
