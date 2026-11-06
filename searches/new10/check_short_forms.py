#!/usr/bin/env python3
"""Check short_forms.txt: each polynomial is APN, quadratic, and has the ortho-derivative invariant of its class."""
import sys, os, csv, hashlib
here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, "..", "..")
sys.path.insert(0, os.path.join(root, "lib"))
from apnsearch.invariants import ortho_invariant, algebraic_degree
N, MOD = 1024, (1 << 10) | 0b1001
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
inv = {r["inv"]: r["name"] for r in csv.DictReader(open(os.path.join(root, "data", "quadratic_classes.tsv")), delimiter="\t")}
want = {"N1": "cc9135ebccd0fd0d", "N2": "dcffd811adb7bb59", "N3": "0a0ea0bfcbf7b32b", "N4": "fc409cf646cabcb1"}
ok = True
for line in open(os.path.join(here, "short_forms.txt")):
    if line.startswith("#"): continue
    name, T, terms = line.rstrip("\n").split("\t")
    terms = [tuple(map(int, t.split(":"))) for t in terms.split()]
    F = [0] * N
    for x in range(N):
        v = 0
        for c, e in terms: v ^= mul(c, pw(x, e))
        F[x] = v
    du = max(max(sum(1 for x in range(N) if F[x ^ a] ^ F[x] == b) for b in set(F[x ^ a] ^ F[x] for x in range(N))) for a in range(1, N))
    h = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
    good = du == 2 and algebraic_degree(F, 10) == 2 and h == want[name] and int(T) == len(terms)
    ok &= good
    print(name, len(terms), "terms, uniformity", du, "invariant", h, inv.get(h, "?"), "OK" if good else "FAIL")
sys.exit(0 if ok else 1)
