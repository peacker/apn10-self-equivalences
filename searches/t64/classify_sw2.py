#!/usr/bin/env python3
"""Name the fingerprints found by sw2 (sw2out/*.out): one representative per fingerprint, full ortho-derivative
invariant (classify.py hash), compared with apndb/n10. Cache in sw2out/keys.tsv; prints newly named keys.
NEW classes: LUT appended to sw2out/new_luts.txt.
usage: classify_sw2.py"""
import glob, os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from switch64 import lut_of, inv_hash
known = {}
for line in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")):
    t = line.rstrip("\n").split("\t")
    if t[0] != "name": known[t[1]] = t[0]
cache = os.path.join(here, "sw2out", "keys.tsv")
done = {}
if os.path.exists(cache):
    for l in open(cache):
        k, h, name = l.rstrip("\n").split("\t")[:3]; done[k] = (h, name)
fo = open(cache, "a")
for f in sorted(glob.glob(os.path.join(here, "sw2out", "s*.out"))):
    for line in open(f):
        p = line.split()
        if len(p) != 11 or not p[0].startswith("K") or p[0] in done: continue
        V = list(map(int, p[1:])); F = lut_of(V); h = inv_hash(F); name = known.get(h, "NEW")
        done[p[0]] = (h, name)
        fo.write("%s\t%s\t%s\t%s\t%s\n" % (p[0], h, name, os.path.basename(f), " ".join(p[1:]))); fo.flush()
        print("%s %s %s (%s)" % (p[0], h, name, os.path.basename(f)), flush=True)
        if name == "NEW":
            with open(os.path.join(here, "sw2out", "new_luts.txt"), "a") as g:
                g.write(h + "\t" + ",".join(map(str, F.tolist())) + "\n")
