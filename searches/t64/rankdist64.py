#!/usr/bin/env python3
"""Rank distributions of V (span of the 10 Hermitian matrices) and C = V^perp for the 16 known type-64 classes
(seeds64.txt). usage: rankdist64.py"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from herm import rank_f4, bits_to_h, pair, span_basis
def dist(B):
    cnt = {}
    for mask in range(1, 1 << len(B)):
        v = 0
        for k in range(len(B)):
            if (mask >> k) & 1: v ^= B[k]
        r = rank_f4(bits_to_h(v)); cnt[r] = cnt.get(r, 0) + 1
    return dict(sorted(cnt.items()))
for line in open(os.path.join(here, "seeds64.txt")):
    name, *codes = line.rstrip("\n").split("\t")
    V = span_basis([int(c) for c in codes[0].split()] if len(codes) == 1 else [int(c) for c in codes])
    # C = V^perp via linear algebra on 25 bits
    rows = [sum(pair(h, 1 << k) << k for k in range(25)) for h in V]
    piv = {}
    for m in rows:
        for b in sorted(piv, reverse=True):
            if (m >> b) & 1: m ^= piv[b]
        if m:
            b = m.bit_length() - 1
            for b2 in list(piv):
                if (piv[b2] >> b) & 1: piv[b2] ^= m
            piv[b] = m
    C = []
    for f in (k for k in range(25) if k not in piv):
        v = 1 << f
        for b, m in piv.items():
            if (m >> f) & 1: v |= 1 << b
        C.append(v)
    print(f"{name}\tdim V {len(V)} {dist(V)}\tdim C {len(C)} {dist(C)}", flush=True)
