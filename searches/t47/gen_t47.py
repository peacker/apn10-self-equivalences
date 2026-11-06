#!/usr/bin/env python3
"""Forced prefix for type 47 (A = B = omega on F_4^3 + identity on F_2^4), input basis from reorder.py
(e0..e3 = Fix(A), then (e4, e5), (e6, e7), (e8, e9) = F_4-coordinates 0, 1, 2 as 1, omega), output basis of
types10.tsv ((e0, e1), (e2, e3), (e4, e5) = F_4-coordinates of Im(B + I), e6..e9 = Fix(B)).
Fixes beta(Fix A, Fix A) to Gold x^3 on F_16 (../force_47_x3.txt) and the intertwiners T_z = beta(z, .) on the
F_4-part to multiplication by theta_z in F_64 = F_4[t]/(t^3 + omega), theta_{e_i} given as F_4-coordinate triples.
usage: gen_t47.py THETA0 THETA1 THETA2 THETA3 > force.txt   (THETA = three F_4 digits 0..3, e.g. 100 010 001 020)"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
def m4(a, b):  # F_4 = {0, 1, w, w^2 = w + 1}, a = bit0 + bit1 * w
    r = 0
    for i in range(2):
        if (b >> i) & 1: r ^= a << i
    if r & 4: r ^= 0b111          # w^2 = w + 1
    return r
def m64(x, y):  # x, y: lists of 3 F_4 coefficients of 1, t, t^2; t^3 = w
    r = [0] * 5
    for i in range(3):
        for j in range(3): r[i + j] ^= m4(x[i], y[j])
    for k in (4, 3):
        r[k - 3] ^= m4(r[k], 2); r[k] = 0
    return r[:3]
th = [[int(c) for c in s] for s in sys.argv[1:5]]
for line in open(os.path.join(here, "..", "..", "engine", "force_47_x3.txt")): print(line.strip())
for i in range(4):
    for k in range(3):
        for s, scal in ((0, 1), (1, 2)):          # input e_{4+2k+s} = scal * t^k
            u = [0, 0, 0]; u[k] = scal
            v = m64(th[i], u); out = 0
            for c in range(3): out |= v[c] << (2 * c)
            j = 4 + 2 * k + s
            print(j * (j - 1) // 2 + i, out)
