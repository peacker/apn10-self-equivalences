#!/usr/bin/env python3
"""Claim: if F on F_2^{2m} (m odd) is quadratic with self-equivalence (A, B), A = B = J_m (m Jordan blocks of
size 2), then every component c.F with c vanishing on Fix(B) is non-bent.  Check on random members of the
type-24 normal form (m = 5) and of its n = 6 analogue (m = 3); then sample APN members at m = 3 and report
their classes (Kim's class should appear if it lies in this type)."""
import random, sys, os, hashlib
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant
from apnsearch import ccz
MODS = {3: 0b1011, 5: 0b100101}
def setup(m):
    MOD = MODS[m]
    def fm(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> m: a ^= MOD
        return r
    gold = [fm(fm(x, x), x) for x in range(1 << m)]
    return fm, gold
def make(m, rng, fm, gold):
    N = 1 << m
    S = [[rng.randrange(N) for _ in range(m)] for _ in range(m)]
    Dq = [[rng.randrange(N) if i < j else 0 for j in range(m)] for i in range(m)]
    def bil(T, x, y):
        r = 0
        for i in range(m):
            if (x >> i) & 1:
                for j in range(m):
                    if (y >> j) & 1: r ^= T[i][j]
        return r
    mu = lambda a, b: gold[a ^ b] ^ gold[a] ^ gold[b]
    L = []
    for x in range(N * N):
        p, q = x & (N - 1), x >> m
        top = gold[p] ^ bil(S, p, p) ^ mu(p, q)
        bot = bil(Dq, p, p) ^ bil(S, q, p) ^ gold[q]
        L.append(top | (bot << m))
    return np.array(L)
def walsh_bent(L, c, n):
    f = np.array([bin(int(v) & c).count("1") & 1 for v in L])
    w = 1 - 2 * f
    h = 1
    while h < len(w):
        w = w.reshape(-1, 2 * h); a, b = w[:, :h].copy(), w[:, h:].copy()
        w[:, :h], w[:, h:] = a + b, a - b; w = w.reshape(-1); h *= 2
    return np.all(np.abs(w) == (1 << (n // 2)))
def is_apn(L, n):
    N = 1 << n; x = np.arange(N)
    for a in range(1, N):
        if len(np.unique(L[x ^ a] ^ L)) != N // 2: return False
    return True
if __name__ == "__main__":
    rng = random.Random(1)
    for m in (3, 5):
        fm, gold = setup(m); n = 2 * m
        bad = 0
        for _ in range(20):
            L = make(m, rng, fm, gold)
            bad += sum(walsh_bent(L, c, n) for c in range(1, 1 << m))   # c on the top (quotient) coordinate
        print(f"m={m}: bent components among the top {2**m-1} (20 random F): {bad}")
    fm, gold = setup(3); seen = {}
    for t in range(200000):
        L = make(3, rng, fm, gold)
        if not is_apn(L, 6): continue
        inv = hashlib.sha1(repr(ortho_invariant(L, 6)).encode()).hexdigest()[:16]
        if inv not in seen:
            st, info = ccz.ccz_to_permutation_dual(L, 6, quadratic=True)
            seen[inv] = (st, info.get("max_dim_NB"))
            print(inv, seen[inv], "after", t, flush=True)
    print("classes:", seen)
