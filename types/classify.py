#!/usr/bin/env python3
"""Classify n = 10 lookup tables (one per line) by the ortho-derivative invariant and name the known classes
(data/quadratic_classes.tsv, a copy of apndb/n10). Prints: invariant, name (or NEW), line count.
usage: classify.py LUTFILE [--ccz]   (--ccz: run the counting criterion + dual test on NEW classes)"""
import collections, hashlib, os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
known = {}
for line in open(os.path.join(here, "..", "data", "quadratic_classes.tsv")):
    t = line.rstrip("\n").split("\t")
    if t[0] != "name":
        known[t[1]] = t[0]
cnt = collections.Counter(); lut = {}
for line in open(sys.argv[1]):
    if "," not in line:
        continue
    F = np.array(list(map(int, line.strip().split(","))))
    inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
    cnt[inv] += 1; lut.setdefault(inv, F)
for inv, c in cnt.most_common():
    name = known.get(inv, "NEW")
    extra = ""
    if name == "NEW" and "--ccz" in sys.argv:
        from apnsearch import ccz
        st, info = ccz.ccz_to_permutation_dual(lut[inv], 10, quadratic=True)
        extra = f"\tM={info.get('max_dim_NB')}\t{st}"
    print(f"{inv}\t{name}\t{c}{extra}")
