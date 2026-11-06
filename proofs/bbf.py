#!/usr/bin/env python3
"""Types excluded by the bound of Baudrin, Beierle, Felke, Leander, Neumann, Perrin, Stennes (DCC 93(8), 2025,
Corollary 3 and Remark 1): if S o A = B o S with (A, B) a non-trivial affine self-equivalence and S has differential
uniformity at most p, then rank(L_A - I) + rank(L_B - I) >= n - 1. For APN functions of F_2^10 this reads
dim Fix(A) + dim Fix(B) <= 11, so every type with dim Fix(A) + dim Fix(B) >= 12 is empty.
Reads ../types/types10.tsv (columns 6, 7: the columns of A and B as integers; A = I or B = I for the kinds 1,p / p,1)
and prints the excluded type ids on one line (the format of bbf_excluded.txt).
usage: bbf.py [-v]"""
import os, sys
n = 10
here = os.path.dirname(os.path.abspath(__file__))


def rank(cols):
    basis = []
    for v in cols:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)


def dim_fix(cols):                                   # dim ker(M + I)
    return n - rank([c ^ (1 << i) for i, c in enumerate(cols)])


excluded = []
for line in open(os.path.join(here, "..", "types", "types10.tsv")):
    t = line.rstrip("\n").split("\t")
    A = list(map(int, t[5].split(","))); B = list(map(int, t[6].split(",")))
    fa, fb = dim_fix(A), dim_fix(B)
    if "-v" in sys.argv:
        print(t[0], t[2], fa, fb, fa + fb, file=sys.stderr)
    if fa + fb >= 12:
        excluded.append(t[0])
print(" ".join(excluded))
