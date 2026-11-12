#!/usr/bin/env python3
"""Type 24, Lemma "regular": the 155 maps L(x) = lambda x^(2^k) (lambda in F_32^*, 0 <= k < 5) of F_32 (polynomial
x^5 + x^2 + 1, as in c/full24g.c) act regularly on the 155 planes (2-dimensional F_2-subspaces) of F_32.

Checks: there are 155 planes; the 155 maps are distinct, F_2-linear and invertible; they send the plane
<e_0, e_1> = {0, 1, x, x + 1} to 155 distinct planes (so every plane, each exactly once); every plane has a trivial
stabilizer. For G = x^3 and G = x^5 (the two Gold functions of the normal form), each L is an automorphism of G:
G o L = N o G with N = G o L o G^(-1) F_2-linear (G is a permutation of F_32), so the group is the same for both.
usage: python3 -I searches/t24/planes_check.py   (exit status 1 on failure)"""
import itertools, sys

POLY, Q = 0b100101, 32


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1; a <<= 1
        if a & Q:
            a ^= POLY
    return r


def pw(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a); e >>= 1
    return r


def linear(f):
    return all(f[x ^ y] == f[x] ^ f[y] for x in range(Q) for y in range(Q))


planes = {frozenset({0, a, b, a ^ b}) for a, b in itertools.combinations(range(1, Q), 2)}
maps = [[mul(l, pw(x, 1 << k)) for x in range(Q)] for l in range(1, Q) for k in range(5)]
ok = len(planes) == 155 and len({tuple(m) for m in maps}) == 155
ok &= all(linear(m) and len(set(m)) == Q for m in maps)
P0 = frozenset({0, 1, 2, 3})
images = {frozenset(m[v] for v in P0) for m in maps}
ok_orbit = images == planes and len(images) == 155
ok_stab = all(sum(frozenset(m[v] for v in P) == P for m in maps) == 1 for P in planes)
ok_aut = True
for e in (3, 5):
    G = [pw(x, e) for x in range(Q)]
    Ginv = {y: x for x, y in enumerate(G)}
    ok_aut &= len(Ginv) == Q and all(linear([G[m[Ginv[y]]] for y in range(Q)]) for m in maps)
print(f"planes {len(planes)}, maps {len(maps)} (distinct, linear, invertible: {ok}); images of <e0, e1>: "
      f"{len(images)} distinct planes = all planes: {ok_orbit}; every stabilizer trivial: {ok_stab}; "
      f"automorphisms of x^3 and x^5 (N linear): {ok_aut}")
ok = ok and ok_orbit and ok_stab and ok_aut
print("REGULAR ACTION OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
