#!/usr/bin/env python3
"""Numerical check of the odd-rank lemma for type 56 (n = 10, A = 5 order-3 blocks, B = 2 blocks + I_6).
Builds the space of alternating bilinear maps beta: V x V -> W with beta(Ax, Ay) = B beta(x, y), samples random
members and checks that some a != 0 has dim ker beta_w(a, .) >= 8, where beta_w is the part of beta in the
B-moving 4-dim subspace W_w = Im(B + I). Then the 6 remaining coordinates cannot cut the kernel down to {0, a}.
usage: check56.py [TYPE] [SAMPLES]   (TYPE row of types10.tsv)"""
import os, sys, random
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
t = sys.argv[1] if len(sys.argv) > 1 else "56"
S = int(sys.argv[2]) if len(sys.argv) > 2 else 200
row = [l.rstrip("\n").split("\t") for l in open(os.path.join(here, "..", "..", "types", "types10.tsv")) if l.split("\t")[0] == t][0]
n = 10
def mat(cols):  # column j = image of e_j, bits = rows
    M = np.zeros((n, n), dtype=np.uint8)
    for j, c in enumerate(map(int, cols.split(","))):
        for i in range(n): M[i, j] = (c >> i) & 1
    return M
A, B = mat(row[5]), mat(row[6])
pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
# unknowns: beta(e_i, e_j) in F_2^10 for i < j -> 450 bits; constraint beta(Ae_i, Ae_j) = B beta(e_i, e_j)
idx = {(i, j, k): (pairs.index((i, j))) * n + k for (i, j) in pairs for k in range(n)}
N = len(pairs) * n
def bvec(u, v):  # coefficient matrix (n x N) of beta(u, v) for vectors u, v
    R = np.zeros((n, N), dtype=np.uint8)
    for i in range(n):
        if not u[i]: continue
        for j in range(n):
            if not v[j] or i == j: continue
            p = (min(i, j), max(i, j))
            for k in range(n): R[k, idx[(p[0], p[1], k)]] ^= 1
    return R
rows = []
for (i, j) in pairs:
    L = bvec(A[:, i], A[:, j])
    Rr = np.zeros((n, N), dtype=np.uint8)
    for k in range(n):
        for l in range(n):
            if B[k, l]: Rr[k, idx[(i, j, l)]] ^= 1
    rows.append((L ^ Rr))
C = np.vstack(rows) % 2
def nullspace(M):
    M = M.copy() % 2; r, c = M.shape; piv = []; row_ = 0
    for col in range(c):
        p = next((x for x in range(row_, r) if M[x, col]), None)
        if p is None: continue
        M[[row_, p]] = M[[p, row_]]
        for x in range(r):
            if x != row_ and M[x, col]: M[x] ^= M[row_]
        piv.append(col); row_ += 1
        if row_ == r: break
    free = [x for x in range(c) if x not in piv]; basis = []
    for f in free:
        v = np.zeros(c, dtype=np.uint8); v[f] = 1
        for k, pc in enumerate(piv):
            if M[k, f]: v[pc] = 1
        basis.append(v)
    return np.array(basis)
def rank(M):
    M = M.copy() % 2; r = 0
    for col in range(M.shape[1]):
        p = next((x for x in range(r, M.shape[0]) if M[x, col]), None)
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for x in range(M.shape[0]):
            if x != r and M[x, col]: M[x] ^= M[r]
        r += 1
    return r
K = nullspace(C)
print(f"type {t}: invariant beta space dim {len(K)}")
# projection to W_w = Im(B + I) along Fix(B): rows of P give coordinates on Im(B+I)
I = np.eye(n, dtype=np.uint8); Nb = (B ^ I) % 2
dimw = rank(Nb)
P = Nb  # beta_w(a, x) = (B + I) beta(a, x); kernel of beta_w(a, .) = kernel of (B+I) beta(a, .)
vec = [np.array([(a >> i) & 1 for i in range(n)], dtype=np.uint8) for a in range(1 << n)]
rng = random.Random(1); worst = []
for s in range(S):
    coef = np.array([rng.randint(0, 1) for _ in range(len(K))], dtype=np.uint8)
    beta = (coef @ K) % 2
    T = np.zeros((n, n, n), dtype=np.uint8)  # T[i, j] = beta(e_i, e_j)
    for (i, j) in pairs:
        for k in range(n): T[i, j, k] = T[j, i, k] = beta[idx[(i, j, k)]]
    best = 0
    for a in range(1, 1 << n):
        u = vec[a]
        M = (np.tensordot(u, T, axes=(0, 0)) % 2).T  # n x n: column j = beta(a, e_j)
        kd = n - rank((P @ M) % 2)
        best = max(best, kd)
        if best >= 8: break
    worst.append(best)
print("max over a of dim ker beta_w(a,.), per sample: min", min(worst), "max", max(worst), f"(W_w dim {dimw})")
