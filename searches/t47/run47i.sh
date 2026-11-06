#!/bin/bash
# run47i.sh MODE CLASS N: type 47 stage 2 with stage2i. MODE p: prefix group N of CLASS (all X representatives
# with the same first three echelon columns); MODE x: the single X slice N (rank-2 classes). Output outi/TAG.out/.done
m=$1; c=$2; n=$3
if [ "$m" = p ]; then tag=$c.p$n; env=PREFIX=$n; else tag=$c.$n; env=XONLY=$n; fi
[ -f outi/$tag.done ] && exit 0
s0=$(date +%s)
env $env LIMIT=250 ./polite.sh ./stage2i cls/$c.txt > outi/$tag.out 2> outi/$tag.err
echo "rc $? secs $(( $(date +%s)-s0 ))" > outi/$tag.done
