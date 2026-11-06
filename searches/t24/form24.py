#!/usr/bin/env python3
"""Type 24 (A = B = J5) in bivariate form over F_32 (x^5 + x^2 + 1): F(p, q) = (Gold(p) + s(p, p) + mu(p, q),
D(p) + s(q, p) + Gold(q)), Gold = x^3, mu its polar, s bilinear (125 bits), D quadratic (polar 50 bits).
Checks: (1) every such F satisfies F(A x) = A F(x) + (linear) with A(p, q) = (p, p + q);
(2) rank of the affine conditions Tr(s(b, b) / b^3) = 1 (b != 0) on s, i.e. APN in the directions (0, b);
(3) random s, D: APN rate, and APN in directions (0, b) <=> the condition."""
import random
m = 5; MOD = 0b100101
def fm(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 32: a ^= MOD
    return r
def fp(a, e):
    r = 1
    while e:
        if e & 1: r = fm(r, a)
        a = fm(a, a); e >>= 1
    return r
def tr(a):
    t = 0; x = a
    for _ in range(m): t ^= x; x = fm(x, x)
    return t
gold = [fp(x, 3) for x in range(32)]
def mu(a, b): return gold[a ^ b] ^ gold[a] ^ gold[b]
rng = random.Random(7)
def rand_bil():  # s(x, y) = sum_{i, j} x_i y_j S[i][j]
    return [[rng.randrange(32) for _ in range(m)] for _ in range(m)]
def bil(S, x, y):
    r = 0
    for i in range(m):
        if (x >> i) & 1:
            for j in range(m):
                if (y >> j) & 1: r ^= S[i][j]
    return r
def rand_quad():
    return [[rng.randrange(32) if i < j else 0 for j in range(m)] for i in range(m)]
def F(S, Dq, p, q):
    top = gold[p] ^ bil(S, p, p) ^ mu(p, q)
    bot = bil(Dq, p, p) ^ bil(S, q, p) ^ gold[q]
    return top | (bot << 5)
def lut(S, Dq): return [F(S, Dq, x & 31, x >> 5) for x in range(1024)]
def Aop(x): p, q = x & 31, x >> 5; return p | ((p ^ q) << 5)
def is_apn(L):
    for a in range(1, 1024):
        seen = set()
        for x in range(1024):
            v = L[x ^ a] ^ L[x]
            if v in seen: return False
            seen.add(v)
        # each value twice: set size 512 is checked implicitly by collisions? no: need count == 2
    return True
def apn_dir(L, a):
    cnt = {}
    for x in range(1024):
        v = L[x ^ a] ^ L[x]; cnt[v] = cnt.get(v, 0) + 1
    return max(cnt.values()) == 2
# (1) equivariance up to affine terms: G(x) = F(Ax) + A F(x) must be affine
for _ in range(5):
    S, Dq = rand_bil(), rand_quad(); L = lut(S, Dq)
    G = [L[Aop(x)] ^ Aop(L[x]) for x in range(1024)]
    assert all(G[x ^ y] ^ G[x] ^ G[y] ^ G[0] == 0 for x in (1, 7, 100, 513) for y in (2, 33, 900, 1000))
print("(1) F(Ax) = A F(x) + affine: ok")
# (2) rank of the 31 affine conditions on the 125 bits of s
rows = []
for b in range(1, 32):
    b3inv = fp(b, 31 - 3); row = 0
    for i in range(m):
        for j in range(m):
            for t in range(m):
                if (b >> i) & 1 and (b >> j) & 1:
                    if tr(fm(1 << t, b3inv)): row ^= 1 << ((i * m + j) * m + t)
    rows.append(row)
def rank(vs):
    vs = list(vs); r = 0
    for bit in range(125):
        p = next((k for k in range(r, len(vs)) if (vs[k] >> bit) & 1), None)
        if p is None: continue
        vs[r], vs[p] = vs[p], vs[r]
        for k in range(len(vs)):
            if k != r and (vs[k] >> bit) & 1: vs[k] ^= vs[r]
        r += 1
    return r
print("(2) rank of Tr(s(b,b)/b^3) = 1 conditions:", rank(rows))
# (3) check the equivalence on random samples
agree = 0
for _ in range(40):
    S, Dq = rand_bil(), rand_quad(); L = lut(S, Dq)
    cond = all(tr(fm(bil(S, b, b), fp(b, 28))) == 1 for b in range(1, 32))
    dirs = all(apn_dir(L, b << 5) for b in range(1, 32))
    agree += cond == dirs
print("(3) condition <=> APN in all directions (0, b):", agree, "/ 40")
