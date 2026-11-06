#!/usr/bin/env python3
"""Enumerate the linear self-equivalence types (A, B) of prime order p (p = 2: involutions, j Jordan blocks of size 2) for quadratic functions of F_2^n, as in
Beierle-Leander (IEEE TIT 2022, Sec. 4.1): types (p, p) with ord A = ord B = p, (p, 1) with B = I, (1, p) with A = I.
An element of order p of GL(n, 2) is, up to conjugacy, a block-diagonal matrix of companion matrices of irreducible
factors of Phi_p (each of degree ord_p(2)), completed by an identity block; it is given by the multiset of factors.
Pairs are reduced under power-similarity (A, B) ~ (A^r, B^r), r coprime to p, which permutes the factors of Phi_p
(the factor with root z goes to the factor with root z^r). A and B are conjugated independently.
Prints one line per type: id, p, kind, factors of A, factors of B, A columns, B columns (as for sym_setup.py).
usage: enumerate.py n [--primes 3,5,...] > types.tsv
"""
import itertools, sys

n = int(sys.argv[1])
primes = None
if "--primes" in sys.argv:
    primes = list(map(int, sys.argv[sys.argv.index("--primes") + 1].split(",")))


def pmul(a, b):                      # GF(2)[x] polynomials as integers
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1; b >>= 1
    return r


def pmod(a, m):
    dm = m.bit_length() - 1
    while a and a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a


def ordp(p):
    k, x = 1, 2 % p
    while x != 1:
        x = x * 2 % p; k += 1
    return k


def factors_phi(p):
    """irreducible factors of Phi_p over GF(2), each labelled by the cyclotomic coset of exponents of its roots."""
    k = ordp(p)
    # work in GF(2^k) via a primitive polynomial found by search; roots of Phi_p are g^((2^k - 1) / p * j)
    K = 1 << k
    for m in range(K | 1, 2 * K, 2):
        # test primitivity of m (degree k): x has order 2^k - 1 modulo m
        x, ok = 1, True
        for e in range(1, K - 1):
            x = pmod(x << 1, m)
            if x == 1:
                ok = False; break
        if ok and pmod(x << 1, m) == 1:
            break
    def fmul(a, b):
        return pmod(pmul(a, b), m)
    def fpow(a, e):
        r = 1
        while e:
            if e & 1:
                r = fmul(r, a)
            a = fmul(a, a); e >>= 1
        return r
    z = fpow(2, (K - 1) // p)          # a primitive p-th root of unity
    seen, facs = set(), []
    for j in range(1, p):
        if j in seen:
            continue
        coset = []
        t = j
        while t not in coset:
            coset.append(t); t = t * 2 % p
        seen.update(coset)
        # minimal polynomial prod (X - z^c), coefficients in GF(2^k) -> must lie in GF(2)
        poly = [1]                      # list of field elements, poly[i] = coeff of X^i
        for c in coset:
            r = fpow(z, c)
            new = [0] * (len(poly) + 1)
            for i, co in enumerate(poly):
                new[i + 1] ^= co
                new[i] ^= fmul(co, r)
            poly = new
        assert all(co in (0, 1) for co in poly)
        facs.append((min(coset), sum(co << i for i, co in enumerate(poly)), frozenset(coset)))
    return k, facs


def companion(poly):
    """columns of the companion matrix of a monic polynomial of degree k (X e_i = e_{i+1}, X e_{k-1} = -sum c_i e_i)"""
    k = poly.bit_length() - 1
    cols = [1 << (i + 1) for i in range(k - 1)]
    cols.append(poly & ((1 << k) - 1))
    return cols


def block_matrix(blocks):
    cols, off = [], 0
    for poly in blocks:
        for c in companion(poly):
            cols.append(c << off)
        off += poly.bit_length() - 1
    while off < n:
        cols.append(1 << off); off += 1
    return cols


out = []
tid = 0
# p = 2: involutions, up to conjugacy j Jordan blocks [[0,1],[1,0]] (companion of (X+1)^2 = X^2 + 1), j = 1..n/2
inv = list(range(1, n // 2 + 1))
for kind, (CA, CB) in (("p,p", (inv, inv)), ("p,1", (inv, [0])), ("1,p", ([0], inv))):
    for ja in CA:
        for jb in CB:
            A = block_matrix([0b101] * ja); B = block_matrix([0b101] * jb)
            out.append((tid, 2, kind, (f"J{ja}",) if ja else (), (f"J{jb}",) if jb else (), A, B)); tid += 1
for p in ([q for q in (primes or []) if q != 2] if primes else
          [q for q in range(3, 1 << n) if all(q % r for r in range(2, q)) and ordp(q) <= n]):
    k, facs = factors_phi(p)
    labels = [f[0] for f in facs]                 # coset leader of each factor
    poly_of = {f[0]: f[1] for f in facs}
    coset_of = {f[0]: f[2] for f in facs}
    def leader(e):
        e %= p
        for l in labels:
            if e in coset_of[l]:
                return l
    classes = [()]                                # multisets (sorted tuples of leaders), total degree <= n
    for r in range(1, n // k + 1):
        classes += [c for c in itertools.combinations_with_replacement(sorted(labels), r)]
    nontriv = [c for c in classes if c]
    def power(c, r):
        return tuple(sorted(leader(l * r) for l in c))
    units = [r for r in range(1, p) if r % p]
    def canon(ca, cb):
        return min((power(ca, r), power(cb, r)) for r in units)
    seen = set()
    for kind, (CA, CB) in (("p,p", (nontriv, nontriv)), ("p,1", (nontriv, [()])), ("1,p", ([()], nontriv))):
        for ca in CA:
            for cb in CB:
                key = canon(ca, cb)
                if (kind, key) in seen:
                    continue
                seen.add((kind, key))
                A = block_matrix([poly_of[l] for l in key[0]])
                B = block_matrix([poly_of[l] for l in key[1]])
                out.append((tid, p, kind, key[0], key[1], A, B)); tid += 1
for t in out:
    print(f"{t[0]}\t{t[1]}\t{t[2]}\t{','.join(map(str, t[3])) or '-'}\t{','.join(map(str, t[4])) or '-'}\t"
          f"{','.join(map(str, t[5]))}\t{','.join(map(str, t[6]))}")
import collections
cnt = collections.Counter((t[1], t[2]) for t in out)
print(f"# {len(out)} types: " + ", ".join(f"p={p} {k}: {c}" for (p, k), c in sorted(cnt.items())), file=sys.stderr)
print(f"# by kind: {dict(collections.Counter(t[2] for t in out))}", file=sys.stderr)
