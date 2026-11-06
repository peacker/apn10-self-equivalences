#!/bin/bash
# run47r.sh CLASS X: type 47 helper on a second server (job list taken in reverse); skips slices finished here
# or claimed by the main server (claimed/TAG, pushed by the local sync script); politely (nice 19, load gate)
c=$1; x=$2; tag=$c.$x; [ -f out2/$tag.done ] && exit 0; [ -f claimed/$tag ] && exit 0
s0=$(date +%s)
XONLY=$x LIMIT=250 ./polite.sh ./stage2g cls/$c.txt > out2/$tag.out 2> out2/$tag.err
echo "rc $? secs $(( $(date +%s)-s0 ))" > out2/$tag.done
