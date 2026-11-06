#!/usr/bin/env python3
"""Name poly64 fingerprints (K line + L line with the LUT) with the full ortho-derivative invariant; cache in
polyout/keys.tsv; NEW classes' LUTs to polyout/new_luts.txt. usage: classify_poly.py FILES..."""
import os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from switch64 import inv_hash
known = {}
for line in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")):
    t = line.rstrip("\n").split("\t")
    if t[0] != "name": known[t[1]] = t[0]
os.makedirs(os.path.join(here, "polyout"), exist_ok=True)
cache = os.path.join(here, "polyout", "keys.tsv"); done = {}
if os.path.exists(cache):
    for l in open(cache): k, h, name = l.rstrip("\n").split("\t")[:3]; done[k] = (h, name)
fo = open(cache, "a")
for f in sys.argv[1:]:
    lines = open(f).read().split("\n")
    for i, line in enumerate(lines):
        if not line.startswith("K"): continue
        k = line.split()[0]
        if k in done or i + 1 >= len(lines) or not lines[i + 1].startswith("L"): continue
        F = np.array(list(map(int, lines[i + 1][1:].split(","))), dtype=np.int64)
        h = inv_hash(F); name = known.get(h, "NEW"); done[k] = (h, name)
        fo.write("%s\t%s\t%s\t%s\n" % (k, h, name, line)); fo.flush()
        print("%s %s %s | %s" % (k, h, name, line), flush=True)
        if name == "NEW":
            with open(os.path.join(here, "polyout", "new_luts.txt"), "a") as g: g.write(h + "\t" + lines[i + 1][1:] + "\n")
