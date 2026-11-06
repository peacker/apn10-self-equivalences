#!/usr/bin/env python3
"""LUTs of the n = 10 functions F_l, F_1, F_2, F_3 of Helleseth-Ichanska-Kaleyski (ePrint 2026/2401, Sec. IV),
for both primitive beta in F_4, modulus x^10 + x^3 + 1."""
N, MOD = 1024, 1024 | 0b1001
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
fams = {"Fl": [(0, 3), (1, 9), (2, 96), (0, 288), (1, 33), (2, 264), (2, 528)],
        "F1": [(0, 3), (1, 9), (2, 96), (0, 288), (1, 66), (1, 132), (0, 264), (0, 528)],
        "F2": [(0, 3), (1, 129), (2, 96), (0, 36), (1, 66), (0, 132), (0, 264), (0, 528)],
        "F3": [(0, 3), (1, 129), (2, 96), (0, 36), (2, 33), (0, 66), (1, 132), (1, 264), (1, 528)]}
for b in (pw(2, 341), pw(2, 682)):
    for name, terms in fams.items():
        lut = []
        for x in range(N):
            v = 0
            for e, k in terms: v ^= mul(pw(b, e), pw(x, k))
            lut.append(v)
        print(",".join(map(str, lut)))
