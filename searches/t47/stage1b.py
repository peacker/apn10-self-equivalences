#!/usr/bin/env python3
"""Type 47 stage 1, part 2: for each anisotropic Hermitian span S (aniso_classes.pkl from stage1.py) compute
Stab(S) in GL_3(F_4) (numpy, all 181440 matrices), its orbits on the F_4-subspaces V of F_4^3 (V = column space of
G, acted on by conj(P)^{-1}), and write the -force file cls/c<S>_<V>.txt with beta on U x U:
  beta(x, y) = G^T (conj x cross conj y) in U' (bits 0..5)  +  sum_k Tr(x^* M_k y) e_{6+k}   (M_k: basis of S)
Checks equivariance beta(wx, wy) = w beta(x, y) and alternation before writing."""
import itertools, os, pickle, sys
import numpy as np
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "stage1.py")).read().split("# anisotropic subspaces")[0])
classes = pickle.load(open("aniso_classes.pkl", "rb"))
MT = np.array(M4, dtype=np.uint8); CJ = np.array([cj(a) for a in range(4)], dtype=np.uint8)
G_all = np.array(GL, dtype=np.uint8)                      # (181440, 3, 3)
def bmul(X, Y):  # batched F_4 matmul, X (m,3,3), Y (m,3,3) or (3,3)
    if Y.ndim == 2: Y = np.broadcast_to(Y, X.shape)
    R = np.zeros(X.shape, dtype=np.uint8)
    for k in range(3): R ^= MT[X[:, :, k][:, :, None], Y[:, k, :][:, None, :]]
    return R
Pstar = CJ[np.transpose(G_all, (0, 2, 1))]
def code(M):  # Hermitian matrix (m,3,3) -> 9-bit index as in stage1.herm
    c = M[:, 0, 0].astype(np.int32) | (M[:, 1, 1].astype(np.int32) << 1) | (M[:, 2, 2].astype(np.int32) << 2)
    c |= (M[:, 0, 1].astype(np.int32) << 3) | (M[:, 0, 2].astype(np.int32) << 5) | (M[:, 1, 2].astype(np.int32) << 7)
    return c
assert all(code(np.array([herm[i]], dtype=np.uint8))[0] == i for i in range(512))
# F_4-subspaces of F_4^3: by their set of vectors (as 6-bit ints, coordinate k = bits 2k, 2k+1)
vecs = list(itertools.product(range(4), repeat=3))
def vint(v): return v[0] | (v[1] << 2) | (v[2] << 4)
def fspan(basis):
    s = set()
    for co in itertools.product(range(4), repeat=len(basis)):
        v = [0, 0, 0]
        for c, b in zip(co, basis):
            for i in range(3): v[i] ^= M4[c][b[i]]
        s.add(tuple(v))
    return frozenset(s)
subs = {frozenset({(0, 0, 0)})}
for v in vecs:
    if any(v): subs.add(fspan([v]))
for v, w in itertools.combinations([v for v in vecs if any(v)], 2):
    sp = fspan([v, w])
    if len(sp) == 16: subs.add(sp)
subs.add(frozenset(vecs)); subs = list(subs)
print("subspaces of F_4^3:", len(subs), file=sys.stderr)
def basis_of(sp):
    b = []
    for v in sorted(sp):
        if any(v) and len(fspan(b + [v])) > len(fspan(b)) if b else any(v): b.append(v)
    return b
def cross(a, b):
    return [M4[a[1]][b[2]] ^ M4[a[2]][b[1]], M4[a[2]][b[0]] ^ M4[a[0]][b[2]], M4[a[0]][b[1]] ^ M4[a[1]][b[0]]]
U = [(b & 3, (b >> 2) & 3, (b >> 4) & 3) for b in range(64)]   # bit vector of U -> F_4^3 (bits 2k, 2k+1 = coord k)
def beta_val(G, Ms, x, y):
    xc = [cj(t) for t in x]; yc = [cj(t) for t in y]; cr = cross(xc, yc)
    out = 0
    for k in range(3):   # (G^T cr)_k = sum_i G[i][k] cr[i]
        s = 0
        for i in range(3): s ^= M4[G[i][k]][cr[i]]
        out |= s << (2 * k)
    for k, Mi in enumerate(Ms):
        M = herm[Mi]; s = 0
        for i in range(3):
            for j in range(3): s ^= M4[M4[cj(x[i])][M[i][j]]][y[j]]
        out |= (s >> 1) << (6 + k)          # Tr(a) = bit 1 of a (Tr(1) = 0, Tr(w) = Tr(w^2) = 1)
    return out
W2 = lambda x: tuple(M4[2][t] for t in x)
def omega_out(v):
    r = 0
    for k in range(3): r |= M4[2][(v >> (2 * k)) & 3] << (2 * k)
    return r | (v & 0b1111000000)
os.makedirs("cls", exist_ok=True)
summary = []
for ci, (dim, osz, rows) in enumerate(classes):
    S = span(rows); Sarr = np.zeros(512, dtype=bool); Sarr[list(S)] = True
    ok = np.ones(len(GL), dtype=bool)
    for r in rows:
        R = bmul(bmul(Pstar, np.array(herm[r], dtype=np.uint8)), G_all)
        ok &= Sarr[code(R)]
    stab = G_all[ok]
    assert len(stab) * osz == len(GL), (len(stab), osz)
    # orbits of {conj(P)} on subspaces
    sid = {s: i for i, s in enumerate(subs)}; seenv = set()
    for s in sorted(subs, key=lambda s: (len(s), sorted(s))):
        if sid[s] in seenv: continue
        orb = set()
        for P in stab:
            Pc = [[cj(int(P[i][j])) for j in range(3)] for i in range(3)]
            img = frozenset(tuple(_dot(Pc[i], v) for i in range(3)) for v in s)
            orb.add(sid[img])
        seenv |= orb
        Vb = basis_of(s) if len(s) > 1 else []
        G = [[0] * 3 for _ in range(3)]
        for k, v in enumerate(Vb):
            for i in range(3): G[i][k] = v[i]
        Ms = list(rows) + [0] * (4 - len(rows))
        # checks: alternating and equivariant on U x U
        for xb in range(64):
            assert beta_val(G, Ms, U[xb], U[xb]) == 0
        for xb, yb in [(1, 2), (5, 9), (17, 40), (63, 33), (12, 48)]:
            assert beta_val(G, Ms, W2(U[xb]), W2(U[yb])) == omega_out(beta_val(G, Ms, U[xb], U[yb]))
        name = f"cls/c{ci:02d}_{len(Vb)}_{len(seenv)}.txt"
        with open(name, "w") as f:
            for j in range(6):
                for i in range(j):
                    f.write(f"{j * (j - 1) // 2 + i} {beta_val(G, Ms, U[1 << i], U[1 << j])}\n")
        summary.append((name, dim, osz, len(stab), len(Vb), len(orb)))
for s in summary: print(*s)
print(len(summary), "classes", file=sys.stderr)
