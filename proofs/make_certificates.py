#!/usr/bin/env python3
"""Write proofs/certificates10.jsonl: one exclusion certificate per type excluded by proof (280 types).

Each type is assigned the first criterion that applies, in the order of Table III of the paper:
  W0    W(A, B) = 0
  C     compression (Proposition "compression") with D a sum of the invariant subspaces ker g(A)^j, j = 1, 2
  R     restriction lemma: D invariant (k = dim D >= 3) with dim E_W(D) = k - 1
  Q     order three lemma: 2 dim Fix(B) < dim U, U = ker(A^2 + A + I)
  WONG  compression with D from the second Wong sequence (wong10.py; same random seed)
  ODD   odd dimension over F_4: A^2 + A + I = 0 (A = omega on F_4^5) and B != I
The witnesses are found with the code of proofs10.py and wong10.py (same search order), and recorded.

One JSON object per line:
  type, p, criterion, A, B   (A, B: lists of 10 column bitmasks, column i = A e_i, as in types/types10.tsv)
  W      basis of W(A, B) as in W/w_ID.txt: per basis vector the 45 values beta(e_i, e_j) (10-bit integers),
         pairs ordered (0,1), (0,2), (1,2), (0,3), ... (by j, then i). Bit r of the value of pair q is unknown
         number 10 q + r of the linear system (450 unknowns).
  a, D   (C, WONG) the vector a and a basis of D (bitmasks); "a_in_D"; "dim_E" = dim span{w_k(a, d_j)}
  D      (R) a basis of D; "dim_E" = dim span{w_k(d, e)}
  dim_U, dim_FixB  (Q)
  invariant  true when D is claimed A-invariant (C, R)
The checker proofs/check_certificates.py verifies every certificate without this code.
usage: python3 -I proofs/make_certificates.py [WDIR]   (default WDIR = W, written by types/make_w.sh)
"""
import json, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402
from proofs10 import (rows, load_w, apply_vec, images, span_basis, in_span, invariant_family, kernel,  # noqa: E402
                      poly_cols, n, N)
import wong10  # noqa: E402

ORDER = ("W0", "C", "R", "Q", "WONG", "ODD")


def read_w_values(path):
    """Basis of W as stored in w_ID.txt: list of 45-value lists."""
    L = open(path).read().split()
    d = int(L[1]); P = n * (n - 1) // 2; pos = 2; out = []
    for _ in range(d):
        pos += 2
        out.append([int(v) for v in L[pos:pos + P]]); pos += P
    return out


def find_C(T, fam):
    for a in range(1, N):
        M = apply_vec(T, a)
        for D in fam:
            S = span_basis(images(M, D))
            need = len(D) - (1 if in_span(a, D) else 0)
            if len(S) < need:
                return a, D, len(S)
    return None


def find_R(T, fam):
    for D in fam:
        k = len(D)
        if k < 3:
            continue
        vals = [images(apply_vec(T, x), D) for x in D]
        E = span_basis(np.unique(np.concatenate(vals)))
        if len(E) == k - 1:
            return D, len(E)
    return None


def find_Q(p, A, B):
    if p != 3:
        return None
    U = kernel(poly_cols(0b111, A)); FB = kernel(poly_cols(0b11, B))
    if U and 2 * len(FB) < len(U):
        return len(U), len(FB)
    return None


def find_wong(T):
    rng = random.Random(1)                      # as in wong10.py
    for a in range(1, N):
        gap, D, S = wong10.search(T, a, 4, rng)
        if gap > 0:
            return a, D, len(S)
    return None


def is_odd(A, B):
    A2 = [wong10.mat_apply(A, c) for c in A]
    return all(A2[i] ^ A[i] ^ (1 << i) == 0 for i in range(n)) and B != [1 << i for i in range(n)]


def main():
    wdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "W")
    out = []
    for tid, r in rows.items():
        A = list(map(int, r[5].split(","))); B = list(map(int, r[6].split(",")))
        p = int(r[1])
        path = os.path.join(wdir, "w_%s.txt" % tid)
        T = load_w(path)
        Wv = read_w_values(path)
        cert = {"type": int(tid), "p": p, "criterion": None, "A": A, "B": B, "W": Wv}
        fam = invariant_family(A)
        if T.shape[0] == 0:
            cert["criterion"] = "W0"
        elif (c := find_C(T, fam)):
            a, D, s = c
            cert.update(criterion="C", a=a, D=list(map(int, D)), a_in_D=in_span(a, list(D)), dim_E=s, invariant=True)
        elif (c := find_R(T, fam)):
            D, e = c
            cert.update(criterion="R", D=list(map(int, D)), dim_E=e, invariant=True)
        elif (c := find_Q(p, A, B)):
            cert.update(criterion="Q", dim_U=c[0], dim_FixB=c[1])
        elif (c := find_wong(T)):
            a, D, s = c
            cert.update(criterion="WONG", a=a, D=list(map(int, D)), a_in_D=in_span(a, list(D)), dim_E=s,
                        invariant=False)
        elif is_odd(A, B):
            cert.update(criterion="ODD")
        else:
            print("type %s: no proof criterion" % tid, file=sys.stderr, flush=True)
            continue
        out.append(cert)
        print("type %s: %s" % (tid, cert["criterion"]), file=sys.stderr, flush=True)
    with open(os.path.join(HERE, "certificates10.jsonl"), "w") as f:
        for c in out:
            f.write(json.dumps(c, separators=(",", ":")) + "\n")
    cnt = {k: sum(c["criterion"] == k for c in out) for k in ORDER}
    print("certificates: %d  %s" % (len(out), " ".join("%s %d" % kv for kv in cnt.items())))


if __name__ == "__main__":
    main()
