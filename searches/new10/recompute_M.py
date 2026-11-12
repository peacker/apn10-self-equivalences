#!/usr/bin/env python3
"""Recompute M (the largest dimension of a linear space of linear forms c, all giving non-bent components c.F,
together with 0) exhaustively, for the four new classes.

Sources: the univariate short forms of short_forms.txt (field x^10 + x^3 + 1, F_4 = {0, 1, 236, 237}) and the
lookup tables of data/quadratic_classes.tsv (column M is compared). A component c.F of a quadratic F is non-bent iff
the alternating form c.beta (beta the polar form) is degenerate. The largest subspace of N = {non-bent c} + {0} is
found by a depth-first search over increasing generators c_1 < c_2 < ... with c_{k+1} not in span(c_1..c_k): every
subspace has such a basis (take c_{k+1} = the smallest element outside the span), so the search is exhaustive.
Expected: N1 3, N2 4, N3 4, N4 4 (as in the paper, Table "new classes").
usage: python3 -I searches/new10/recompute_M.py [--all]   (--all: also the 25 published classes of the tsv)
Exit status 1 on a mismatch."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MOD = 0x409
EXPECTED = {"N1": 3, "N2": 4, "N3": 4, "N4": 4}


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1; a <<= 1
        if a & 0x400:
            a ^= MOD
    return r


def pw(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a); e >>= 1
    return r


def rank2(rows):
    rows = list(rows); r = 0
    for bit in range(10):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def max_M(lut):
    beta = lambda x, y: lut[x ^ y] ^ lut[x] ^ lut[y] ^ lut[0]
    B = [[beta(1 << i, 1 << j) for j in range(10)] for i in range(10)]
    nonbent = set()
    for c in range(1, 1024):
        rows = [sum(((bin(B[i][j] & c).count("1") & 1) << j) for j in range(10)) for i in range(10)]
        if rank2(rows) < 10:
            nonbent.add(c)
    cand = sorted(nonbent)
    best = [0]

    def ext(span, dim, last):
        best[0] = max(best[0], dim)
        for c in cand:
            if c <= last or c in span:
                continue
            new = span | {s ^ c for s in span}
            if all(v == 0 or v in nonbent for v in new):
                ext(new, dim + 1, c)

    ext({0}, 0, 0)
    return len(nonbent), best[0]


assert mul(236, 236) ^ 236 ^ 1 == 0
ok = True
for line in open(os.path.join(HERE, "short_forms.txt")):
    if line.startswith("#") or not line.strip():
        continue
    name, nt, terms = line.rstrip("\n").split("\t")
    T = [tuple(map(int, t.split(":"))) for t in terms.split()]
    lut = [0] * 1024
    for x in range(1024):
        v = 0
        for c, e in T:
            v ^= mul(c, pw(x, e))
        lut[x] = v
    nb, M = max_M(lut)
    good = M == EXPECTED[name]
    ok &= good
    print(f"short form {name}: non-bent components {nb}, M = {M} (expected {EXPECTED[name]}){'' if good else '  MISMATCH'}")

rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(ROOT, "data", "quadratic_classes.tsv"))]
hdr = rows[0]; iM = hdr.index("M"); iL = hdr.index("lut")
k = 0
for r in rows[1:]:
    is_new = r[0].startswith("New")
    if not is_new and "--all" not in sys.argv:
        continue
    lut = list(map(int, r[iL].split(",")))
    nb, M = max_M(lut)
    good = str(M) == r[iM]
    if is_new:
        k += 1
        good &= M == EXPECTED[f"N{k}"]
    ok &= good
    print(f"tsv {r[0][:60]}: non-bent components {nb}, M = {M} (tsv column M = {r[iM]}){'' if good else '  MISMATCH'}")
print("M OK" if ok else "M MISMATCH")
sys.exit(0 if ok else 1)
