#!/usr/bin/env python3
"""Rewrite A in an input basis where C(A)-invariant subspaces come first: for an involution (N = A + I) a basis of
Im N, extended to a basis of Fix(A) = ker N, then a complement; for odd order a basis of Fix(A), then the rest.
Prints the columns of A' = P^{-1} A P (for sym_setup.py --amat), the dimensions of the invariant subspaces (the
levels where a checkpoint can use the whole C(A)), and checks A' A' ... has the same order.
usage: reorder.py TYPE_ID [TYPES_FILE [n]]"""
import os, sys
tid = sys.argv[1]
path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "types10.tsv")
row = [l.rstrip("\n").split("\t") for l in open(path) if l.split("\t")[0] == tid][0]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 10
A = list(map(int, row[5].split(",")))
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
def extend(basis, pool):
    for v in pool:
        if rank(basis + [v]) > len(basis): basis = basis + [v]
    return basis
N_ = [A[i] ^ (1 << i) for i in range(n)]                  # N = A + I (columns)
fix = [v for v in range(1, 1 << n) if app(A, v) == v]     # Fix(A) = ker N
img = [app(N_, v) for v in range(1, 1 << n)]
img = [v for v in img if v]
B1 = extend([], sorted(set(img))) if int(row[1]) == 2 else []
B2 = extend(B1, fix)
P = extend(B2, [1 << i for i in range(n)])
assert len(P) == n
# coordinates in basis P: solve by brute force over the span (n = 10: 1024 combinations)
coord = {}
for t in range(1 << n):
    v = 0
    for i in range(n):
        if (t >> i) & 1: v ^= P[i]
    coord[v] = t
Ap = [coord[app(A, P[i])] for i in range(n)]
print(",".join(map(str, Ap)))
print(f"# type {tid}: dim Im N = {len(B1)}, dim Fix(A) = {len(B2)} (checkpoint levels with the whole C(A))", file=sys.stderr)
