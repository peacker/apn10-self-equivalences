#!/bin/bash
# run47n.sh CLASS X: type 47 stage 2 with the dedicated solver (stage2g.c) for one class and one X slice, politely
c=$1; x=$2; tag=$c.$x; [ -f out2/$tag.done ] && exit 0
s0=$(date +%s)
XONLY=$x LIMIT=250 ./polite.sh ./stage2g cls/$c.txt > out2/$tag.out 2> out2/$tag.err
echo "rc $? secs $(( $(date +%s)-s0 ))" > out2/$tag.done
