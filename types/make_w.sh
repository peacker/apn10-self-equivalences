#!/bin/bash
# make_w.sh [OUTDIR] [TYPES...]: write OUTDIR/w_ID.txt (default OUTDIR = ../W, all 299 types) with sym_setup.py:
# the space W of equivariant alternating bilinear maps for type ID, followed by the orbit representatives of
# beta(e_0, e_1), a basis of the centraliser of B and a basis of the centraliser of A (sym10 input format).
# sym_setup.py is deterministic (fixed seed), so the files equal the ones used by the searches (engine/W3).
here="$(cd "$(dirname "$0")" && pwd)"; out="${1:-$here/../W}"; shift 2>/dev/null
mkdir -p "$out"
ids="$*"; [ -z "$ids" ] && ids=$(cut -f1 "$here/types10.tsv")
for id in $ids; do
  line=$(awk -F'\t' -v i="$id" '$1==i' "$here/types10.tsv")
  A=$(echo "$line" | cut -f6); B=$(echo "$line" | cut -f7)
  python3 "$here/sym_setup.py" 10 --amat "$A" --bmat "$B" -o "$out/w_$id.txt" > /dev/null
done
