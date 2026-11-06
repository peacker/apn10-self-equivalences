#!/usr/bin/env python3
"""Forced prefix for sym10: beta on Fix(A) x Fix(A) set to the bilinear form of a Gold function x^(2^k + 1) on
F_{2^f}, f = dim Fix(A). Assumes the input basis from reorder.py (the first f basis vectors span Fix(A)); the
output values are written in a basis of Fix(B) (B from types10.tsv, canonical coordinates). Valid when
C(A) and C(B) act on Fix(A), Fix(B) as the full GL(f) (e.g. A = B = J5, or semisimple A, B with equal fixed
dimension), so that every APN restriction can be moved to the Gold representative.
usage: gen_force.py TYPE_ID K MODULUS > force.txt   (pair index = j (j - 1) / 2 + i for i < j < f)"""
import os, sys
tid, k, mod = sys.argv[1], int(sys.argv[2]), int(sys.argv[3], 0)
row = [l.rstrip("\n").split("\t") for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "types", "types10.tsv")) if l.split("\t")[0] == tid][0]
n = 10; B = list(map(int, row[6].split(",")))
def app(M, v):
    r = 0
    for i in range(n):
        if (v >> i) & 1: r ^= M[i]
    return r
def rank(vs):
    rows = list(vs); r = 0
    for b in range(n):
        p = next((i for i in range(r, len(rows)) if (rows[i] >> b) & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> b) & 1: rows[i] ^= rows[r]
        r += 1
    return r
fixB = []
for v in range(1, 1 << n):
    if app(B, v) == v and rank(fixB + [v]) > len(fixB): fixB.append(v)
f = mod.bit_length() - 1                         # field degree
assert len(fixB) == f, (len(fixB), f)
F = 1 << f
def fmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & F: a ^= mod
    return r
def fpow(a, e):
    r = 1
    while e:
        if e & 1: r = fmul(r, a)
        a = fmul(a, a); e >>= 1
    return r
q = 1 << k
for j in range(f):
    for i in range(j):
        x, y = 1 << i, 1 << j                    # polynomial basis 1, g, g^2, ... of F_{2^f}
        val = fmul(fpow(x, q), y) ^ fmul(x, fpow(y, q))       # bilinear form of x^(q+1)
        out = 0
        for t in range(f):
            if (val >> t) & 1: out ^= fixB[t]
        print(j * (j - 1) // 2 + i, out)
