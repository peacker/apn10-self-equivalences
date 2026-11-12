#!/usr/bin/env python3
"""Dimension six check for Section "Type 24" of the paper: a x^5 + a x^33 + x^40 with a = g^54 in F_8 is APN and
EA-equivalent to Kim's function x^3 + x^10 + g x^24 on F_64, g a root of x^6 + x^4 + x^3 + x + 1 (primitive).
The equivalence (L1, L2) is computed by lib/apnsearch/equivalence.ea_equivalence and re-checked by
check_certificate. Both functions are also compared by the differential spectrum of the ortho-derivative.
usage: python3 -I searches/t24/kim6_f8.py   (exit status 1 on failure)"""
import os, sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "lib"))
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402

MOD = 0b1011011                  # x^6 + x^4 + x^3 + x + 1
N = 64


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1; a <<= 1
        if a & N:
            a ^= MOD
    return r


def pw(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a); e >>= 1
    return r


def lut(terms):                  # sum of c x^e
    return [sum_xor(mul(c, pw(x, e)) for c, e in terms) for x in range(N)]


def sum_xor(it):
    r = 0
    for v in it:
        r ^= v
    return r


def is_apn(L):
    return all(max(Counter(L[x] ^ L[x ^ a] for x in range(N)).values()) <= 2 for a in range(1, N))


def ortho(L):                    # ortho-derivative of a quadratic APN function
    pi = [0] * N
    for a in range(1, N):
        img = {L[x] ^ L[x ^ a] ^ L[a] ^ L[0] for x in range(N)}
        pi[a] = next(c for c in range(1, N) if all(bin(c & v).count("1") % 2 == 0 for v in img))
    return pi


def dspec(T):
    C = Counter()
    for a in range(1, N):
        cnt = Counter(T[x] ^ T[x ^ a] for x in range(N))
        for b in range(N):
            C[cnt.get(b, 0)] += 1
    return tuple(sorted(C.items()))


g = 2
ok = True
prim = len({pw(g, i) for i in range(63)}) == 63
a = pw(g, 54)
in_f8 = pw(a, 8) == a
kim = lut([(1, 3), (1, 10), (g, 24)])
F8f = lut([(a, 5), (a, 33), (1, 40)])
print("g primitive:", prim, "| a = g^54 =", a, "in F_8:", in_f8)
print("Kim x^3 + x^10 + g x^24 APN:", is_apn(kim))
print("a x^5 + a x^33 + x^40 APN:", is_apn(F8f))
same = dspec(ortho(kim)) == dspec(ortho(F8f))
print("ortho-derivative differential spectra equal:", same)
r = ea_equivalence(kim, F8f, 6)
cert = r is not None and check_certificate(kim, F8f, 6, r[0], r[1])
print("EA-equivalent (explicit L1, L2, certificate checked):", cert)
ok = prim and in_f8 and is_apn(kim) and is_apn(F8f) and same and cert
print("OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
