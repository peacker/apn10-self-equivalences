#!/bin/bash
# run_all.sh: phase 1 = 67 search jobs (4 types x 16 shards + 170, 171, 231), phase 2 = EA chunks; then FPGAP_DONE.
cd "$(dirname "$0")"; mkdir -p out ea
{ for t in 129 175 176 180; do for s in $(seq 0 15); do echo "$t $s"; done; done; echo 170; echo 171; echo 231; } > jobs1.txt
xargs -P 67 -L 1 ./search_job.sh < jobs1.txt
{ for f in out/*.done; do tag=$(basename $f .done)
    case $tag in t_176.*) k=12;; t_175.*|t_180.*) k=12;; t_170|t_171) k=24;; *) k=1;; esac
    for c in $(seq 0 $((k-1))); do echo "$tag $c $k"; done; done; } > jobs2.txt
xargs -P 180 -L 1 ./ea_job.sh < jobs2.txt
touch FPGAP_DONE
