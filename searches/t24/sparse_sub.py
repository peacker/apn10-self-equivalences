#!/usr/bin/env python3
"""Sparse quadratic polynomials on F_{2^{2m}} with coefficients in the subfield F_{2^m}: sum a_e x^e over
e = 2^i + 2^j (0 <= i < j < 2m), at most K nonzero terms, first coefficient 1. Prints APN classes found
(ortho invariant, counting M, permutation test).  usage: sparse_sub.py m K"""
import itertools, sys, os, hashlib
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant
from apnsearch import ccz, fast
m, K = int(sys.argv[1]), int(sys.argv[2]); n = 2 * m; N = 1 << n
POLY = {6: 0b1000011, 10: 0b10000001001}[n]
exp = [0] * (2 * N); log = [0] * N; x = 1
for k in range(N - 1):
    exp[k] = x; log[x] = k; x <<= 1
    if x >> n: x ^= POLY
for k in range(N - 1, 2 * N): exp[k] = exp[k - (N - 1)]
exp = np.array(exp); log = np.array(log)
X = np.arange(N)
def mono(e):  # x^e as a table
    out = np.zeros(N, dtype=np.int64); nz = X != 0
    out[nz] = exp[(log[X[nz]] * e) % (N - 1)]; return out
def scal(a, T):
    if a == 0: return np.zeros(N, dtype=np.int64)
    out = np.zeros(N, dtype=np.int64); nz = T != 0
    out[nz] = exp[(log[a] + log[T[nz]]) % (N - 1)]; return out
sub = [y for y in range(N) if y and exp[(log[y] * (1 << m)) % (N - 1)] == y]   # F_{2^m}^*
E = sorted({(1 << i) + (1 << j) for i in range(n) for j in range(i + 1, n)})
# up to x -> x^2 (cyclotomic shift) we may take the first exponent as a class representative
reps = sorted({min(((e << s) | (e >> (n - s))) & (N - 1) for s in range(n)) for e in E})
M = {e: mono(e) for e in E}
seen = {}; tested = 0
for k in range(1, K + 1):
    for e0 in reps:
        for rest in itertools.combinations([e for e in E if e > e0], k - 1):
            for coefs in itertools.product(sub, repeat=k - 1):
                T = M[e0].copy()
                for a, e in zip(coefs, rest): T ^= scal(a, M[e])
                tested += 1
                if not fast.is_apn(T): continue
                inv = hashlib.sha1(repr(ortho_invariant(T, n)).encode()).hexdigest()[:16]
                if inv in seen: continue
                st, info = ccz.ccz_to_permutation_dual(T, n, quadratic=True)
                seen[inv] = (st, info.get("max_dim_NB"), [e0] + list(rest), [1] + [int(a) for a in coefs])
                print(inv, seen[inv], flush=True)
print("tested", tested, "classes", len(seen))
