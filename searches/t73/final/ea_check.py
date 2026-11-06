#!/usr/bin/env python3
"""EA-equivalence check of the type-73 leaves: group all printed leaves (out5/*.out, K lines) by the
ortho-derivative invariant; inside each group test every leaf against the first one, and test the first one
against the apndb representative with the same invariant (certificates re-checked).
usage: ea_check.py [group_index]"""
import collections, glob, hashlib, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402
known = {}
for line in open(os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")):
    t = line.rstrip("\n").split("\t")
    if t[0] != "name":
        known[t[1]] = (t[0], t)
hdr = open(os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")).readline().rstrip("\n").split("\t")
groups = collections.OrderedDict()
for f in sorted(glob.glob(os.path.join(here, "out5", "*.out"))):
    for line in open(f):
        if not line.startswith("K"): continue
        F = np.array(list(map(int, line.split()[1].split(","))))
        inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
        groups.setdefault(inv, []).append((os.path.basename(f), F))
print(f"{sum(len(v) for v in groups.values())} printed leaves, {len(groups)} invariant groups", flush=True)
lutcol = next((i for i, h in enumerate(hdr) if "lut" in h.lower() or "table" in h.lower()), None)
sel = [int(a) for a in sys.argv[1:]] or range(len(groups))
for gi, (inv, items) in enumerate(groups.items()):
    if gi not in sel: continue
    name = known.get(inv, ("?", None))[0]
    t0 = time.time(); R = items[0][1]; bad = 0
    for src, F in items[1:]:
        c = ea_equivalence(R, F, 10)
        if c is None or not check_certificate(R, F, 10, *c): bad += 1; print(f"  NOT EQUIVALENT: {src}", flush=True)
    ref = "no apndb LUT"
    if lutcol is not None and inv in known:
        D = np.array(list(map(int, known[inv][1][lutcol].split(","))))
        c = ea_equivalence(D, R, 10)
        ref = "equivalent to apndb representative" if (c is not None and check_certificate(D, R, 10, *c)) else "NOT equivalent to apndb representative"
    print(f"group {gi} {inv} {name}: {len(items)} leaves, {len(items) - 1 - bad} equivalent to the first, {bad} not; {ref}; {time.time() - t0:.0f} s", flush=True)
