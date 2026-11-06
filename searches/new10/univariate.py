#!/usr/bin/env python3
"""Conjugate the type-73 matrix A (= B) to the Frobenius x -> x^4 of F_2^10 (modulus x^10 + x^3 + 1), transform each
LUT G = P F P^-1 (then G(x^4) = G(x)^4) and print its univariate polynomial.
usage: univariate.py LUTFILE TYPES_FILE"""
import sys
n, N, MOD = 10, 1024, (1 << 10) | 0b1001
row = [l.rstrip("\n").split("\t") for l in open(sys.argv[2]) if l.split("\t")[0] == "73"][0]
A = list(map(int, row[5].split(","))); assert row[5] == row[6]
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & N: a ^= MOD
    return r
def pw(a, e):
    r = 1
    while e:
        if e & 1: r = mul(r, a)
        a = mul(a, a); e >>= 1
    return r
Fr = [pw(1 << i, 4) for i in range(n)]
def app(M, v):
    r = 0
    for i in range(n):
        if (v >> i) & 1: r ^= M[i]
    return r
def rank(vs):
    rows = list(vs); r = 0
    for b in range(n):
        p = next((i for i in range(r, len(rows)) if (rows[i] >> b) & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> b) & 1: rows[i] ^= rows[r]
        r += 1
    return r
def adapted(M):
    fix = [v for v in range(1, N) if app(M, v) == v]
    phi = [v for v in range(1, N) if (lambda w: w[0] ^ w[1] ^ w[2] ^ w[3] ^ w[4])(
        [v, app(M, v), app(M, app(M, v)), app(M, app(M, app(M, v))), app(M, app(M, app(M, app(M, v))))]) == 0]
    bas = []
    for v in fix:
        if rank(bas + [v]) > len(bas): bas.append(v)
    assert len(bas) == 2
    for v in phi:
        orb = [v, app(M, v), app(M, app(M, v)), app(M, app(M, app(M, v)))]
        if rank(bas + orb) == len(bas) + 4: bas += orb
        if len(bas) == 10: break
    assert len(bas) == 10
    return bas
bA, bF = adapted(A), adapted(Fr)
# P: v -> w with P(bA[i]) = bF[i]
coordA = {}
for t in range(N):
    v = 0
    for i in range(n):
        if (t >> i) & 1: v ^= bA[i]
    coordA[v] = t
def P(v):
    t = coordA[v]; w = 0
    for i in range(n):
        if (t >> i) & 1: w ^= bF[i]
    return w
Pinv = {P(v): v for v in range(N)}
assert all(P(app(A, v)) == app(Fr, P(v)) for v in range(N))
for line in open(sys.argv[1]):
    F = list(map(int, line.split(",")))
    G = [P(F[Pinv[x]]) for x in range(N)]
    terms = []
    for k in range(1, N - 1):
        c = 0
        for x in range(1, N):
            if G[x]: c ^= mul(G[x], pw(x, N - 1 - k))
        if c and bin(k).count("1") == 2: terms.append((k, c))   # quadratic part only (EA-invariant here)
    c0 = G[0]
    print(f"quadratic terms (exponent: coefficient, coefficient in F_4 = {{0,1,{pw(2,341)},{pw(2,682)}}}?):")
    print("  " + ", ".join(f"x^{k}: {c}{'*' if pw(c, 4) == c else ''}" for k, c in terms))
