#!/usr/bin/env python3
"""Type 47 stage 1: classes of the U x U part of beta (U = F_4^3, canonical basis e0..e5, coordinate k = bits 2k, 2k+1
as 1, omega). beta(x, y) for x, y in U = (C(x, y), H-part) with
  C(x, y) = G^T (conj x  cross  conj y)  in U' = F_4^3   (G in M_3(F_4); the output group leaves only V = col G),
  h_k(x, y) = Tr(x^* M_k y), M_k Hermitian           in Z' = F_2^4 (the output group leaves only S = span M_k).
Input P in GL_3(F_4) acts by S -> P^* S P, V -> conj(P)^{-1} V. Necessary: S anisotropic (x^* M x = 0 for all M in S
only for x = 0). Prints one -force file per class (S, V) as cls_<id>.txt and a summary line per class."""
import itertools, sys
M4 = [[0] * 4 for _ in range(4)]
for a in range(4):
    for b in range(4):
        r = 0
        for i in range(2):
            if (b >> i) & 1: r ^= a << i
        if r & 4: r ^= 7
        M4[a][b] = r
def cj(a): return M4[a][a]
INV = {a: b for a in range(1, 4) for b in range(1, 4) if M4[a][b] == 1}
def mm(X, Y): return tuple(tuple(_dot(X[i], [Y[k][j] for k in range(3)]) for j in range(3)) for i in range(3))
def _dot(u, v):
    r = 0
    for a, b in zip(u, v): r ^= M4[a][b]
    return r
def ct(X): return tuple(tuple(cj(X[j][i]) for j in range(3)) for i in range(3))
def cjm(X): return tuple(tuple(cj(x) for x in row) for row in X)
def det(X):
    a = X
    t = lambda i, j, k: M4[M4[a[0][i]][a[1][j]]][a[2][k]]
    return t(0,1,2) ^ t(0,2,1) ^ t(1,0,2) ^ t(1,2,0) ^ t(2,0,1) ^ t(2,1,0)
# Hermitian matrices: 9 bits
herm = []
for bits in range(512):
    d = [(bits >> i) & 1 for i in range(3)]
    o = [(bits >> (3 + 2 * i)) & 3 for i in range(3)]
    X = [[0] * 3 for _ in range(3)]
    for i in range(3): X[i][i] = d[i]
    for t, (i, j) in enumerate(((0, 1), (0, 2), (1, 2))): X[i][j] = o[t]; X[j][i] = cj(o[t])
    herm.append(tuple(map(tuple, X)))
hidx = {h: i for i, h in enumerate(herm)}
pts = []  # 21 projective points of PG(2, 4): first nonzero coordinate 1
for v in itertools.product(range(4), repeat=3):
    if any(v) and v[next(i for i in range(3) if v[i])] == 1: pts.append(v)
def qf(Mi, x):  # x^* M x
    M = herm[Mi]; r = 0
    for i in range(3):
        for j in range(3): r ^= M4[M4[cj(x[i])][M[i][j]]][x[j]]
    return r
zeroset = [sum(1 << p for p, x in enumerate(pts) if qf(Mi, x) == 0) for Mi in range(512)]
ALL = (1 << 21) - 1
def span(basis):
    s = {0}
    for b in basis: s |= {x ^ b for x in s}
    return frozenset(s)
# group: all of GL_3(F_4)
mats = [m for m in itertools.product(range(4), repeat=9)]
GL = []
for m in mats:
    X = (m[0:3], m[3:6], m[6:9])
    if det(X): GL.append(X)
print(f"|GL_3(F_4)| = {len(GL)}", file=sys.stderr)
def act_h(P, Mi): return hidx[mm(mm(ct(P), herm[Mi]), P)]
# anisotropic subspaces of dim <= 4 (each subspace once, by reduced echelon form), then orbits under generators
import random
rng = random.Random(0)
gens = [rng.choice(GL) for _ in range(4)]
hperm = [[act_h(P, i) for i in range(512)] for P in gens]
def rref_spaces(dim, nb=9):
    for piv in itertools.combinations(range(nb), dim):        # pivot = leading (highest) bit of each row
        free = []
        for r, p in enumerate(piv):
            free.append([b for b in range(p) if b not in piv])
        nfree = sum(len(f) for f in free)
        for code in range(1 << nfree):
            rows = []; c = code
            for r, p in enumerate(piv):
                v = 1 << p
                for b in free[r]:
                    if c & 1: v |= 1 << b
                    c >>= 1
                rows.append(v)
            yield rows
classes = []
for dim in range(1, 5):
    aniso = {}
    for rows in rref_spaces(dim):
        sp = span(rows); z = ALL
        for m in sp:
            if m: z &= zeroset[m]
            if not z: break
        if z == 0: aniso[frozenset(sp)] = rows
    keys = list(aniso); idx = {k: i for i, k in enumerate(keys)}; par = list(range(len(keys)))
    def f(i):
        while par[i] != i: par[i] = par[par[i]]; i = par[i]
        return i
    for k in keys:
        for hp in hperm:
            t = frozenset(hp[x] for x in k); a, b = f(idx[k]), f(idx[t])
            if a != b: par[a] = b
    orbs = {}
    for k in keys: orbs.setdefault(f(idx[k]), []).append(k)
    print(f"dim {dim}: anisotropic subspaces {len(keys)}, orbits {len(orbs)}, sizes {sorted(len(v) for v in orbs.values())}", file=sys.stderr)
    for o in orbs.values():
        rep = min(o, key=lambda k: sorted(k)); classes.append((dim, len(o), aniso[rep]))
import pickle
pickle.dump(classes, open("aniso_classes.pkl", "wb"))
