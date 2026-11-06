#!/usr/bin/env python3
"""Shrunk-subspace (compression) search with the second Wong sequence, for one direction a at a time.
M_a = {beta_t(a, .) : t in W}. For M0 in M_a: W_0 = 0, W_{i+1} = M_a(M0^{-1}(W_i)) (images under every basis
element of M_a), D = M0^{-1}(W*). Then M_a(D) is contained in W* (when M0 has maximal rank), and the type is empty
if dim D - [a in D] > dim span M_a(D) (checked directly, so any output is a proof whatever M0 was).
usage: wong10.py WDIR TYPES..."""
import os, random, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from proofs10 import load_w, apply_vec, images, span_basis, in_span, n, N

def mat_apply(cols, v):
    r = 0
    for j in range(n):
        if (v >> j) & 1: r ^= cols[j]
    return r

def preimage(cols, Wb):          # {v : M0 v in span Wb}
    # solve via enumeration over basis: build linear map v -> (M0 v mod span Wb)
    red = []
    for c in cols:
        x = c
        for b in sorted(Wb, reverse=True):
            x = min(x, x ^ b)
        red.append(x)
    # kernel of v -> sum v_j red[j]
    piv = []; ker = []
    for j in range(n):
        img, pre = red[j], 1 << j
        for bi, bp in piv:
            if img ^ bi < img: img ^= bi; pre ^= bp
        if img: piv.append((img, pre)); piv.sort(reverse=True)
        else: ker.append(pre)
    return span_basis(ker)

def search(T, a, tries=20, rng=random):
    M = apply_vec(T, a)                    # (d, n): column j of each basis element
    d = M.shape[0]
    best = None
    for _ in range(tries):
        t = np.array([rng.getrandbits(1) for _ in range(d)], dtype=bool)
        M0 = [int(x) for x in (np.bitwise_xor.reduce(M[t], axis=0) if t.any() else np.zeros(n, dtype=np.uint16))]
        Wb = []
        while True:
            D = preimage(M0, Wb)
            S = span_basis(images(M, D)) if D else []
            if len(S) == len(Wb): break
            Wb = S
        D = preimage(M0, Wb); S = span_basis(images(M, D)) if D else []
        gap = len(D) - (1 if D and in_span(a, D) else 0) - len(S)
        if best is None or gap > best[0]: best = (gap, D, S)
    return best

if __name__ == "__main__":
    wdir = sys.argv[1]
    for tid in sys.argv[2:]:
        T = load_w("%s/w_%s.txt" % (wdir, tid)); rng = random.Random(1)
        res = None
        for a in range(1, N):
            gap, D, S = search(T, a, 4, rng)
            if gap > 0: res = (a, D, S); break
        if res:
            a, D, S = res
            print("%s\tproof: a=%d, dim D=%d (a %s D), dim beta(a, D)=%d; D=%s" % (tid, a, len(D), "in" if in_span(a, D) else "not in", len(S), D), flush=True)
        else:
            print("%s\tno compression found" % tid, flush=True)
