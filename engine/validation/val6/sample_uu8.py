#!/usr/bin/env python3
"""Type 47, canonical basis (U = e0..e5 as F_4^3, Z = e6..e9): draw random members of W, keep the U x U part
(pairs i < j < 6) when q(u) = beta(u, A u) is nonzero for all u in U \\ 0 (the Hermitian part is anisotropic, a
necessary condition), and print it as a -force file.  usage: sample_uu.py WFILE SEED > force.txt"""
import sys, random
lines = open(sys.argv[1]).read().split("\n")
n, d = map(int, lines[0].split()); basis = [list(map(int, l.split()[2:])) for l in lines[1:1 + d]]
rng = random.Random(int(sys.argv[2]))
P = n * (n - 1) // 2
pidx = {}
for j in range(n):
    for i in range(j): pidx[(i, j)] = j * (j - 1) // 2 + i
A = [2, 3, 8, 12, 16, 32, 64, 128]
def app(M, v):
    r = 0
    for i in range(n):
        if (v >> i) & 1: r ^= M[i]
    return r
while True:
    v = [0] * P
    for b in basis:
        if rng.random() < 0.5: v = [x ^ y for x, y in zip(v, b)]
    def beta(x, y):
        r = 0
        for i in range(n):
            if not (x >> i) & 1: continue
            for j in range(n):
                if i != j and (y >> j) & 1: r ^= v[pidx[(min(i, j), max(i, j))]]
        return r
    if all(beta(u, app(A, u)) for u in range(1, 16)):
        for j in range(4):
            for i in range(j): print(pidx[(i, j)], v[pidx[(i, j)]])
        break
