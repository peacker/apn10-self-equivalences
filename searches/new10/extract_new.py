#!/usr/bin/env python3
"""For LUTs whose ortho-derivative fingerprint is in the given list: print one LUT per fingerprint, then check
differential uniformity (direct), algebraic degree, M and the dual test, and print the ortho-derivative spectra.
usage: extract_new.py LUTFILE FP1 [FP2 ...]"""
import hashlib, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant, algebraic_degree  # noqa: E402
from apnsearch import ccz  # noqa: E402
want = set(sys.argv[2:]); got = {}
for line in open(sys.argv[1]):
    if "," not in line: continue
    F = np.array(list(map(int, line.strip().split(","))))
    inv = ortho_invariant(F, 10)
    fp = hashlib.sha1(repr(inv).encode()).hexdigest()[:16]
    if fp in want and fp not in got: got[fp] = (F, inv)
for fp, (F, inv) in got.items():
    X = np.arange(1024)
    du = max(np.bincount(F[X ^ a] ^ F, minlength=1024).max() for a in range(1, 1024))
    st, info = ccz.ccz_to_permutation_dual(F, 10, quadratic=True)
    print(f"# {fp}: uniformity {du}, degree {algebraic_degree(F, 10)}, M={info.get('max_dim_NB')}, perm test {st}")
    print(f"# diff {dict(inv[0])}")
    print(",".join(map(str, F)))
