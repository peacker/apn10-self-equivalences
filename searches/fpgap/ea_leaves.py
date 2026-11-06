#!/usr/bin/env python3
"""EA check of printed leaves: ea_leaves.py OUT.gz CHUNK NCHUNKS > RESULT.
Takes every CHUNK-th (mod NCHUNKS) leaf line "K<fp> <lut>" of OUT.gz. For each leaf: ortho-derivative invariant;
the apndb representative with the same invariant; exhaustive EA test (kernels=True, transitive=True for Gold
classes), certificate re-checked. One line per leaf: <sha1(lut)[:16]> <inv> <class>: <equivalent|NOT EQUIVALENT|
NOT IN APNDB> <secs>. A failing leaf's LUT is printed after it on a line starting with "L"."""
import gzip, hashlib, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402
known = {t[1]: t for t in (l.rstrip("\n").split("\t") for l in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")))
         if t[0] != "name"}
path, ci, nc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
reps = {}; k = -1
for line in gzip.open(path, "rt"):
    if not line.startswith("K"): continue
    k += 1
    if k % nc != ci: continue
    s = line.split()[1]; h = hashlib.sha1(s.encode()).hexdigest()[:16]
    F = np.array(list(map(int, s.split(","))))
    t0 = time.time()
    inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
    if inv not in known: res = "NOT IN APNDB"; name = "?"
    else:
        name = known[inv][0]
        D = reps.setdefault(inv, np.array(list(map(int, known[inv][-1].split(",")))))
        c = ea_equivalence(D, F, 10, transitive=name.startswith("Gold"), kernels=True)
        res = "equivalent" if (c is not None and check_certificate(D, F, 10, *c)) else "NOT EQUIVALENT"
    print(f"{h} {inv} {name}: {res} {time.time() - t0:.1f}", flush=True)
    if res != "equivalent": print("L" + s, flush=True)
print("END", k + 1, flush=True)
