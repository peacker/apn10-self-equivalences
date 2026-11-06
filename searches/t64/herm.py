#!/usr/bin/env python3
"""Type 64 (F(omega x) = F(x) + linear on F_2^10 = F_4^5) in Hermitian form.
Each component c.F of such F has polar form Tr_{F4/F2}(x^* H_c y) with H_c in Herm_5(F_4). F is APN iff the
trace dual C = V^perp (V = span H_c, pairing <H, G> = tr(HG) in F_2) has no nonzero element of rank <= 2.
This module: F_4 / GF(2^10) arithmetic, Hermitian matrices of a quadratic map given by its polar form, ranks.
usage: herm.py   (self-test on x^3, x^9, x^3 + Tr(x^9) and a non-APN control)"""
import itertools, sys

MOD = 0x409; n = 10
def gmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> n: a ^= MOD
    return r
def gpow(a, e):
    r = 1
    while e:
        if e & 1: r = gmul(r, a)
        a = gmul(a, a); e >>= 1
    return r
OM = gpow(2, 341)                      # omega in GF(2^10)
assert OM != 1 and gmul(OM, gmul(OM, OM)) == 1
# F_4 = {0, 1, w, w^2} as 0, 1, 2, 3; multiplication table
F4M = [[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]
CONJ = [0, 1, 3, 2]
def f4_to_g(c): return [0, 1, OM, gmul(OM, OM)][c]
# F_4-basis of GF(2^10): greedy
def f4span(vs):
    s = {0}
    for v in vs:
        s = {x ^ gmul(f4_to_g(c), v) for x in s for c in range(4)}
    return s
basis = []
for v in range(1, 1 << n):
    if v not in f4span(basis): basis.append(v)
    if len(basis) == 5: break
def vec(coords):                        # F_4^5 -> GF(2^10)
    r = 0
    for c, b in zip(coords, basis): r ^= gmul(f4_to_g(c), b)
    return r

def herm_of_polar(Bfun):
    """Bfun(x, y) -> 0/1, F_2-bilinear, alternating, with Bfun(w x, w y) = Bfun(x, y). Returns H (5x5 over F_4) with
    Bfun(x, y) = Tr(x^* H y), Tr(z) = z + z^2 on F_4 (Tr(1) = 0, Tr(w) = Tr(w^2) = 1)."""
    H = [[0] * 5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            # h = H_ij; Tr(e_i^* H (lam e_j)) = Tr(h lam): lam = w gives Tr(w h), lam = w^2 gives Tr(w^2 h)
            t1 = Bfun(vec([1 if k == i else 0 for k in range(5)]), vec([2 if k == j else 0 for k in range(5)]))
            t2 = Bfun(vec([1 if k == i else 0 for k in range(5)]), vec([3 if k == j else 0 for k in range(5)]))
            # find h in F_4 with Tr(w h) = t1, Tr(w^2 h) = t2
            for h in range(4):
                tr = lambda z: 0 if z in (0, 1) else 1
                if tr(F4M[2][h]) == t1 and tr(F4M[3][h]) == t2: H[i][j] = h; break
    return H

def rank_f4(M):
    M = [row[:] for row in M]; r = 0
    for c in range(5):
        p = next((i for i in range(r, 5) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = [0, 1, 3, 2][M[r][c]]
        M[r] = [F4M[inv][x] for x in M[r]]
        for i in range(5):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [M[i][k] ^ F4M[f][M[r][k]] for k in range(5)]
        r += 1
    return r

# coordinates of Herm_5(F_4) over F_2 (25 bits): diagonal (5 bits, in F_2) + upper triangle (10 entries, 2 bits)
def h_to_bits(H):
    bits = []
    for i in range(5): bits.append(H[i][i] & 1)
    for i in range(5):
        for j in range(i + 1, 5): bits += [H[i][j] & 1, H[i][j] >> 1]
    return int("".join(map(str, reversed(bits))), 2)
def bits_to_h(v):
    H = [[0] * 5 for _ in range(5)]; k = 0
    for i in range(5): H[i][i] = (v >> k) & 1; k += 1
    for i in range(5):
        for j in range(i + 1, 5):
            h = ((v >> k) & 1) | (((v >> (k + 1)) & 1) << 1); k += 2
            H[i][j] = h; H[j][i] = CONJ[h]
    return H
def pair(u, v):                          # tr(H G) in F_2
    H, G = bits_to_h(u), bits_to_h(v); s = 0
    for i in range(5):
        for k in range(5): s ^= F4M[H[i][k]][G[k][i]]
    assert s in (0, 1)
    return s

def span_basis(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return b

def code_of(F):
    """F: GF(2^10) -> GF(2^10) quadratic with F(w x) = F(x) + linear. Returns (V basis, C = V^perp basis)."""
    Hs = []
    for c in range(n):
        Bf = lambda x, y, c=c: ((F(x ^ y) ^ F(x) ^ F(y) ^ F(0)) >> c) & 1
        Hs.append(h_to_bits(herm_of_polar(Bf)))
    V = span_basis(Hs)
    # C = {G : pair(H, G) = 0 for H in V}: solve over the 25-dim space
    rows = [[pair(h, 1 << k) for k in range(25)] for h in V]
    # null space of rows (len(V) x 25)
    M = [sum(r[k] << k for k in range(25)) for r in rows]
    piv = []; A = []
    for m in M:
        for (pb, pm) in piv:
            if (m >> pb) & 1: m ^= pm
        if m:
            pb = m.bit_length() - 1
            piv = [(b, x ^ m if (x >> pb) & 1 else x) for b, x in piv]
            piv.append((pb, m))
    pivs = {b: m for b, m in piv}
    free = [k for k in range(25) if k not in pivs]
    C = []
    for f in free:
        v = 1 << f
        for b, m in pivs.items():
            if (m >> f) & 1: v |= 1 << b
        C.append(v)
    return V, C

def min_rank(C):
    best = 9; cnt = {}
    for mask in range(1, 1 << len(C)):
        v = 0
        for k in range(len(C)):
            if (mask >> k) & 1: v ^= C[k]
        r = rank_f4(bits_to_h(v)); cnt[r] = cnt.get(r, 0) + 1
        best = min(best, r)
    return best, cnt

def is_apn(F):
    for a in range(1, 1 << n):
        seen = set()
        for x in range(1 << n):
            y = F(x ^ a) ^ F(x)
            if y in seen and False: pass
            seen.add(y)
        if len(seen) != (1 << (n - 1)): return False
    return True

if __name__ == "__main__":
    def tr(z):
        s = 0; x = z
        for _ in range(n): s ^= x; x = gmul(x, x)
        return s
    tests = {"x^3": lambda x: gpow(x, 3), "x^9": lambda x: gpow(x, 9),
             "x^3+Tr(x^9)": lambda x: gpow(x, 3) ^ tr(gpow(x, 9)),
             "x^3+x^9 (control)": lambda x: gpow(x, 3) ^ gpow(x, 9)}
    for name, F in tests.items():
        inv = all(F(gmul(OM, x)) ^ F(x) ^ F(gmul(OM, 0)) == F(gmul(OM, x)) ^ F(x) for x in [0])
        V, C = code_of(F)
        mr, cnt = min_rank(C)
        print("%-20s APN %s  dim V %d  dim C %d  min rank %d  rank counts %s" % (name, is_apn(F), len(V), len(C), mr, dict(sorted(cnt.items()))), flush=True)
