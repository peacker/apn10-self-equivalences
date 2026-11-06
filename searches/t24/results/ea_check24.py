#!/usr/bin/env python3
"""EA check of the type-24 leaves (results/out/*.lut after unpacking t24_out.tar.gz, one leaf per line: tag then LUT, space separated): each leaf
against the apndb representative with the same ortho-derivative invariant (certificate re-checked)."""
import glob, hashlib, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402
tsv = os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")
hdr = open(tsv).readline().rstrip("\n").split("\t")
lutcol = next(i for i, h in enumerate(hdr) if "lut" in h.lower() or "table" in h.lower())
known = {}
for line in open(tsv):
    t = line.rstrip("\n").split("\t")
    if t[0] != "name": known[t[1]] = t
for f in sorted(glob.glob(os.path.join(here, "out", "*.lut"))):
    for line in open(f):
        w = line.split()
        if len(w) < 2: continue
        F = np.array(list(map(int, w[1:] if len(w) > 2 else w[1].split(","))))
        assert len(F) == 1024
        inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
        t0 = time.time()
        if inv not in known: print(f"{os.path.basename(f)} {inv}: NOT IN APNDB", flush=True); continue
        D = np.array(list(map(int, known[inv][lutcol].split(","))))
        c = ea_equivalence(D, F, 10)
        ok = c is not None and check_certificate(D, F, 10, *c)
        print(f"{os.path.basename(f)} {inv} {known[inv][0]}: {'equivalent' if ok else 'NOT equivalent'} to the apndb representative; {time.time() - t0:.0f} s", flush=True)
