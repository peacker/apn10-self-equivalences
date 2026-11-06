#!/usr/bin/env python3
"""Which known n = 10 quadratic APN classes lie in type 64? F has a self-equivalence (A, I) with A of order 3 without
fixed points (up to EA-equivalence) iff its polar form beta has an isometry A of that kind: beta(Ax, Ay) = beta(x, y).
Finds the full group {A : beta(Ax, Ay) = beta(x, y)} by backtracking over the images of the basis vectors.
usage: iso64.py [CLASSES_TSV]"""
import sys
n = 10; N = 1 << n
path = sys.argv[1] if len(sys.argv) > 1 else __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..", "data", "quadratic_classes.tsv")

def isometries(F):
    beta = lambda x, y: F[x ^ y] ^ F[x] ^ F[y] ^ F[0]
    B = [[beta(1 << i, 1 << j) for j in range(n)] for i in range(n)]
    # beta(y, z) for all y, z is costly; precompute row tables lazily
    cache = {}
    def row(y):
        r = cache.get(y)
        if r is None:
            r = [F[y ^ z] ^ F[y] ^ F[z] ^ F[0] for z in range(N)]; cache[y] = r
        return r
    sols = []
    img = []
    def span_has(v, vs):
        b = []
        for x in vs + [v]:
            for c in b: x = min(x, x ^ c)
            if x: b.append(x)
            elif x == 0 and len(b) < len(vs) + 1 and x is not None: pass
        return len(b) < len(vs) + 1
    def rec(k):
        if k == n:
            sols.append(list(img)); return
        rows = [row(img[i]) for i in range(k)]
        for y in range(1, N):
            if all(rows[i][y] == B[i][k] for i in range(k)) and not span_has(y, img):
                img.append(y); rec(k + 1); img.pop()
    rec(0)
    return sols

def app(M, v):
    r = 0
    for i in range(n):
        if (v >> i) & 1: r ^= M[i]
    return r

def order_and_fix(M):
    P = list(M); o = 1
    while P != [1 << i for i in range(n)]:
        P = [app(M, v) for v in P]; o += 1
    fix = sum(1 for v in range(1, N) if app(M, v) == v)
    return o, fix

for line in open(path):
    t = line.rstrip("\n").split("\t")
    if t[0] == "name": continue
    F = list(map(int, t[4].split(",")))
    G = isometries(F)
    kinds = {}
    for M in G:
        o, fx = order_and_fix(M); kinds[(o, fx)] = kinds.get((o, fx), 0) + 1
    in64 = any(o == 3 and fx == 0 for (o, fx) in kinds)
    print("%-28s |Isom|=%-4d orders(order,fixed pts)=%s  TYPE64=%s" % (t[0], len(G), dict(sorted(kinds.items())), "YES" if in64 else "no"), flush=True)
