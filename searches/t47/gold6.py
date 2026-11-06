#!/usr/bin/env python3
"""U x U part of x^3 on F_64 in the basis adapted to A = B = Frobenius x -> x^4 (order 3): input/output bits
0..3 = b1, A b1, b2, A b2 (U = ker(A^2 + A + 1) = F_4^2, coordinate k = bits 2k, 2k+1 as 1, omega), bits 4, 5 =
1, g^21 (Fix(A) = F_4).  Writes the -force style file for stage2g (-DKU=2 -DMZ=2) on stdout.
usage: gold6.py [e]  (exponent, default 3)"""
import sys
e = int(sys.argv[1]) if len(sys.argv) > 1 else 3
MOD = 0b1011011  # x^6 + x^4 + x^3 + x + 1
def mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & 64: a ^= MOD
    return r
def pw(a, k):
    r = 1
    for _ in range(k): r = mul(r, a)
    return r
g = 2
A = lambda x: pw(x, 4)
U = [x for x in range(64) if pw(x, 16) ^ pw(x, 4) ^ x == 0]
assert len(U) == 16
def rank(vs):
    b = []
    for v in vs:
        for w in b: v = min(v, v ^ w)
        if v: b.append(v)
    return len(b)
b1 = U[1]; b2 = next(x for x in U if rank([b1, A(b1), x]) == 3)
basis = [b1, A(b1), b2, A(b2), 1, pw(g, 21)]
assert rank(basis) == 6
# coordinates in the basis
coordof = {}
for m in range(64):
    v = 0
    for i in range(6):
        if m >> i & 1: v ^= basis[i]
    coordof[v] = m
F = lambda x: pw(x, e)
beta = lambda x, y: F(x ^ y) ^ F(x) ^ F(y) ^ F(0)
for j in range(4):
    for i in range(j):
        print(j * (j - 1) // 2 + i, coordof[beta(basis[i], basis[j])])
