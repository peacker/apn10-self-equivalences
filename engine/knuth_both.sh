#!/bin/bash
# knuth_both.sh ID: regenerate W3/w_ID.txt (with Z(B) and Z(A)) and estimate the n = 10 tree with both reductions
id=$1
line=$(awk -F'\t' -v i=$id '$1==i' ../types/types10.tsv); A=$(echo "$line" | cut -f6); B=$(echo "$line" | cut -f7)
nice -n 19 python3 ../types/sym_setup.py 10 --amat "$A" --bmat "$B" -o W3/w_$id.txt > /dev/null
k=$(timeout 900 ./polite.sh ./sym10cp -knuth 2000 -seed $id -symdepth 45 -cp 3 10 < W3/w_$id.txt 2>&1 | grep "^total")
echo -e "$id\t$(echo "$k" | awk '{print $3}' | tr -d ,)\t$(echo "$k" | awk '{print $NF}')" >> knuth_both.tsv
