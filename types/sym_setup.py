#!/usr/bin/env python3
"""Space of quadratic functions of F_2^n with a prescribed linear self-equivalence F(A x) = B F(x) (up to linear
terms), as input for sym10.

A quadratic F (no affine part) is given by its alternating bilinear map beta(x, y) = F(x + y) + F(x) + F(y), i.e.
by the vectors beta(e_i, e_j) in F_2^n for i < j: n * C(n, 2) bits. The self-equivalence means
beta(A x, A y) = B beta(x, y), a homogeneous linear system on these bits. Its solution space W is written in reduced
row echelon form with respect to the search order of sym10 (pairs (i, j) ordered by j, then i; bits inside a pair
in increasing order), so that the value of a pair depends only on basis vectors whose pivot is at or before it.

A and B are given either as field multiplications in GF(2^n) (polynomial basis, modulus MOD):
  --mult ZA ZB        A x = ZA * x, B y = ZB * y   (ZA, ZB field elements as integers)
  --mulexp E K        A = mult by g^E, B = mult by g^(E*K) (g = the generator 2 of the field, i.e. the element x)
or as matrices (--amat / --bmat: n comma-separated column integers, column i = A e_i).
Output: n d, then d lines "pivot_pair pivot_bit v_0 .. v_{P-1}" (P = C(n, 2) pair values of the basis vector), then
a line "r rep_1 .. rep_r": the allowed values of beta(e_0, e_1), one per orbit of the centraliser of B (output
normalisation; the orbits come from 60 random invertible centraliser elements, so they can only be too fine,
never too coarse: the normalisation is always valid), then a line "z" and z lines with the columns of a basis of
the centraliser algebra {X : XB = BX} (used by sym10 -symdepth), then the same for {X : XA = AX} (sym10 -cp).
usage: sym_setup.py n (--mulexp E K | --mult ZA ZB | --amat .. --bmat ..) [--mod 0x409] -o OUT
"""
import argparse, sys

ap = argparse.ArgumentParser()
ap.add_argument("n", type=int)
ap.add_argument("--mod", type=lambda s: int(s, 0), default=0x409)
ap.add_argument("--mult", type=int, nargs=2)
ap.add_argument("--mulexp", type=int, nargs=2)
ap.add_argument("--amat"); ap.add_argument("--bmat")
ap.add_argument("-o", "--out", required=True)
a = ap.parse_args()
n = a.n; N = 1 << n


def fmul(x, y):
    r = 0
    while y:
        if y & 1:
            r ^= x
        y >>= 1
        x <<= 1
        if x & N:
            x ^= a.mod
    return r


def fpow(x, e):
    r = 1
    while e:
        if e & 1:
            r = fmul(r, x)
        x = fmul(x, x); e >>= 1
    return r


if a.mulexp:
    E, K = a.mulexp
    za = fpow(2, E); zb = fpow(2, (E * K) % (N - 1))
    A = [fmul(za, 1 << i) for i in range(n)]
    B = [fmul(zb, 1 << i) for i in range(n)]
elif a.mult:
    A = [fmul(a.mult[0], 1 << i) for i in range(n)]
    B = [fmul(a.mult[1], 1 << i) for i in range(n)]
else:
    A = list(map(int, a.amat.split(","))); B = list(map(int, a.bmat.split(",")))

pairs = [(i, j) for j in range(n) for i in range(j)]            # search order: by j, then i
pidx = {p: k for k, p in enumerate(pairs)}
P = len(pairs)
nv = P * n                                                       # variable (p, r) -> p * n + r


def bit(v, i):
    return (v >> i) & 1


def apply(M, v):                                                 # matrix (columns) times vector
    r = 0
    for i in range(n):
        if bit(v, i):
            r ^= M[i]
    return r


rows = []
for (i, j) in pairs:
    Ai, Aj = A[i], A[j]
    # beta(A e_i, A e_j) = sum_{k<l} (Ai_k Aj_l + Ai_l Aj_k) beta(e_k, e_l)
    coef = []
    for (k, l) in pairs:
        if (bit(Ai, k) & bit(Aj, l)) ^ (bit(Ai, l) & bit(Aj, k)):
            coef.append(pidx[(k, l)])
    for r in range(n):
        eq = 0
        for q in coef:
            eq ^= 1 << (q * n + r)
        # minus B beta(e_i, e_j), bit r: sum_s B[s] bit r * x_{(i,j), s}
        for s in range(n):
            if bit(B[s], r):
                eq ^= 1 << (pidx[(i, j)] * n + s)
        if eq:
            rows.append(eq)

# null space of the system: Gaussian elimination, then free variables
piv = {}
for eq in rows:
    for c, r0 in piv.items():
        if (eq >> c) & 1:
            eq ^= r0
    if not eq:
        continue
    c = eq.bit_length() - 1
    for c2 in list(piv):
        if (piv[c2] >> c) & 1:
            piv[c2] ^= eq
    piv[c] = eq
free = [v for v in range(nv) if v not in piv]
basis = []
for f in free:
    vec = 1 << f
    for c, r0 in piv.items():
        if (r0 >> f) & 1:
            vec |= 1 << c
    basis.append(vec)
# reduced row echelon form with pivot = lowest variable index (search order)
ech = []
for v in basis:
    for e in ech:
        lo = e & -e
        if v & lo:
            v ^= e
    if v:
        lo = v & -v
        ech = [e ^ v if e & lo else e for e in ech]
        ech.append(v)
ech.sort(key=lambda e: (e & -e).bit_length())
d = len(ech)

# output normalisation: invertible maps C commuting with B keep the symmetry (C beta(Ax, Ay) = B C beta(x, y)) and
# APN; so beta(e_0, e_1) can be fixed to one representative per orbit of the centraliser of B on nonzero vectors.
# Orbits: union-find over images under random invertible elements of the centraliser algebra {X : XB = BX}.
import random
rr = random.Random(1)
def centraliser_basis(M):
    """basis of the algebra {X : X M = M X} (X as n columns packed into one integer, column i at bits i*n..)"""
    rows_c = []
    for i in range(n):
        for r in range(n):
            eq = 0
            # (XM) e_i = X (M e_i) = sum_s M[i]_s X e_s ; (MX) e_i = M (X e_i)
            for s2 in range(n):
                if bit(M[i], s2):
                    eq ^= 1 << (s2 * n + r)
            for s2 in range(n):                 # bit r of M (X e_i) = sum_t X[i]_t * bit r of M[t]
                if bit(M[s2], r):
                    eq ^= 1 << (i * n + s2)
            if eq:
                rows_c.append(eq)
    cp = {}
    for eq in rows_c:
        for c, r0 in cp.items():
            if (eq >> c) & 1:
                eq ^= r0
        if not eq:
            continue
        c = eq.bit_length() - 1
        for c2 in list(cp):
            if (cp[c2] >> c) & 1:
                cp[c2] ^= eq
        cp[c] = eq
    basis_c = []
    for f0 in [v for v in range(n * n) if v not in cp]:
        vec = 1 << f0
        for c, r0 in cp.items():
            if (r0 >> f0) & 1:
                vec |= 1 << c
        basis_c.append(vec)
    return basis_c


cbasis = centraliser_basis(B)                   # Z(B): output side
abasis = centraliser_basis(A)                   # Z(A): input side


def cmat_rank(cols):
    rows_ = list(cols); rk = 0
    for b_ in range(n):
        pv_ = next((k for k in range(rk, n) if (rows_[k] >> b_) & 1), None)
        if pv_ is None:
            continue
        rows_[rk], rows_[pv_] = rows_[pv_], rows_[rk]
        for k in range(n):
            if k != rk and (rows_[k] >> b_) & 1:
                rows_[k] ^= rows_[rk]
        rk += 1
    return rk


parent = list(range(N))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
gens = 0
for _ in range(4000):
    if gens >= 60:
        break
    v = 0
    for b_ in cbasis:
        if rr.random() < 0.5:
            v ^= b_
    X = [(v >> (i * n)) & (N - 1) for i in range(n)]
    if cmat_rank(X) < n:
        continue
    gens += 1
    for x in range(1, N):
        y = apply(X, x)
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)
reps = sorted({find(x) for x in range(1, N)})
with open(a.out, "w") as f:
    f.write(f"{n} {d}\n")
    for e in ech:
        pv = (e & -e).bit_length() - 1
        vals = [(e >> (p * n)) & (N - 1) for p in range(P)]
        f.write(f"{pv // n} {pv % n} " + " ".join(map(str, vals)) + "\n")
    f.write(f"{len(reps)} " + " ".join(map(str, reps)) + "\n")
    # basis of the centraliser algebra Z(B) = {X : X B = B X} (for the stabiliser-chain symmetry reduction of
    # sym10): "z" then z lines of n column integers
    f.write(f"{len(cbasis)}\n")
    for v in cbasis:
        f.write(" ".join(str((v >> (i * n)) & (N - 1)) for i in range(n)) + "\n")
    # basis of Z(A) = {X : XA = AX} (input side, checkpoint pruning of sym10): "za" then za lines
    f.write(f"{len(abasis)}\n")
    for v in abasis:
        f.write(" ".join(str((v >> (i * n)) & (N - 1)) for i in range(n)) + "\n")
print(f"n={n}: {len(rows)} equations, d = {d}; beta(e0,e1) orbit reps under C(B): {len(reps)}; pivots per level:",
      [sum(1 for e in ech if pairs[((e & -e).bit_length() - 1) // n][1] == j) for j in range(n)])
