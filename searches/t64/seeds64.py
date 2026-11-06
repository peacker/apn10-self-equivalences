#!/usr/bin/env python3
"""Hermitian form (V, 10 matrices as 25-bit codes) of every known apndb/n10 class that lies in type 64.
For each class: an isometry A of order 3 without fixed points (iso64.py), the F_4-structure omega x = A x, an
F_4-basis b_1..b_5, and H_k from the polar form of component k. Writes seeds64.txt: name, then the 10 codes.
usage: seeds64.py"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from iso64 import isometries, order_and_fix, app
F4M = [[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]
def run():
    out = open(os.path.join(here, "seeds64.txt"), "w")
    for line in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")):
        t = line.rstrip("\n").split("\t")
        if t[0] == "name": continue
        F = list(map(int, t[4].split(",")))
        A = next((M for M in isometries(F) if order_and_fix(M) == (3, 0)), None)
        if A is None: continue
        mul = lambda c, v: [0, v, app(A, v), app(A, v) ^ v][c]      # omega^2 = omega + 1
        def f4span(vs):
            s = {0}
            for v in vs: s = {x ^ mul(c, v) for x in s for c in range(4)}
            return s
        basis = []
        for v in range(1, 1024):
            if v not in f4span(basis): basis.append(v)
            if len(basis) == 5: break
        def vec(coords):
            r = 0
            for c, b in zip(coords, basis): r ^= mul(c, b)
            return r
        e = lambda i, c: vec([c if k == i else 0 for k in range(5)])
        tr = lambda z: 0 if z in (0, 1) else 1
        codes = []
        for k in range(10):
            B = lambda x, y: ((F[x ^ y] ^ F[x] ^ F[y] ^ F[0]) >> k) & 1
            H = [[0] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    t1, t2 = B(e(i, 1), e(j, 2)), B(e(i, 1), e(j, 3))
                    H[i][j] = next(h for h in range(4) if tr(F4M[2][h]) == t1 and tr(F4M[3][h]) == t2)
            bits = [H[i][i] & 1 for i in range(5)]
            for i in range(5):
                for j in range(i + 1, 5): bits += [H[i][j] & 1, H[i][j] >> 1]
            codes.append(sum(b << q for q, b in enumerate(bits)))
        out.write(t[0] + "\t" + " ".join(map(str, codes)) + "\n"); out.flush()
        print(t[0], flush=True)
run()
