#!/usr/bin/env python3
"""EA check of every leaf of the type-73 rerun: distinct LUTs, grouped by ortho-derivative invariant; each distinct
LUT is tested against the apndb representative with the same invariant (certificate re-checked).
Results are appended to ea_all.txt (one line per distinct LUT), so the script can be restarted."""
import glob, hashlib, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402
tsv = os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")
hdr = open(tsv).readline().rstrip("\n").split("\t")
lutcol = next(i for i, h in enumerate(hdr) if "lut" in h.lower() or "table" in h.lower())
known = {t[1]: t for t in (l.rstrip("\n").split("\t") for l in open(tsv)) if t[0] != "name"}
luts = {}
for f in sorted(glob.glob(os.path.join(here, "e?", "out", "a*.out"))):
    for line in open(f):
        if line.startswith("K"): luts.setdefault(line.split()[1], os.path.basename(f))
outp = os.path.join(here, "ea_all.txt")
done = set(l.split()[0] for l in open(outp)) if os.path.exists(outp) else set()
print(f"{len(luts)} distinct LUTs, {len(done)} already checked", flush=True)
reps = {}
with open(outp, "a") as out:
    for s, src in luts.items():
        h = hashlib.sha1(s.encode()).hexdigest()[:16]
        if h in done: continue
        F = np.array(list(map(int, s.split(","))))
        inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
        t0 = time.time()
        if inv not in known: res = "NOT IN APNDB"
        else:
            D = reps.setdefault(inv, np.array(list(map(int, known[inv][lutcol].split(",")))))
            c = ea_equivalence(D, F, 10, transitive=known[inv][0].startswith("Gold"), kernels=True)
            res = "equivalent" if (c is not None and check_certificate(D, F, 10, *c)) else "NOT EQUIVALENT"
        out.write(f"{h} {src} {inv} {known.get(inv, ['?'])[0]}: {res} {time.time() - t0:.0f}s\n"); out.flush()
