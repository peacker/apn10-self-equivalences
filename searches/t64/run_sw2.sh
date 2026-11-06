#!/bin/bash
# run_sw2.sh I: codimension-2 switching from seed line I of seeds64.txt -> sw2out/s$I.out
i=$1; V=$(sed -n "${i}p" seeds64.txt | cut -f2)
./sw2 $V > sw2out/s$i.out 2> sw2out/s$i.err
echo "done $(tail -1 sw2out/s$i.err)" >> sw2out/s$i.err
