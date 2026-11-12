#!/usr/bin/env python3
"""Independent checker of the exclusion certificates proofs/certificates10.jsonl (written by make_certificates.py).

It imports nothing from the repository: the GF(2) linear algebra is reimplemented here on Python integers
(a vector of F_2^k is an int, bit i = coordinate i). For every certificate:

 (i)   from A and B alone, build the homogeneous linear system beta(A e_i, A e_j) = B beta(e_i, e_j), i < j, on
       the 450 unknowns of an alternating bilinear beta (unknown 10 q + r = bit r of beta(e_i, e_j), q the index
       of the pair (i, j) in the order (0,1), (0,2), (1,2), (0,3), ...), and compute its rank r. Check that the
       stored basis of W has 450 - r linearly independent vectors that all solve the system: it spans the full
       solution space W(A, B). Also check that A and B are the matrices of the type in types/types10.tsv,
       that A^p = B^p = I and (A, B) != (I, I).
 (ii)  C, WONG: a != 0, D linearly independent, and dim D - [a in D] > dim span{w_k(a, d_j)} (w_k basis of W,
       d_j basis of D): Proposition "compression". For C, D is A-invariant.
       R: D linearly independent and A-invariant, k = dim D >= 3, dim span{w_k(d_i, d_j)} <= k - 1
       (restriction lemma: then dim E = k - 1 for an APN function, which the lemma forbids).
 (iii) Q: p = 3, dim U and dim Fix(B) recomputed (U = ker(A^2 + A + I)), 2 dim Fix(B) < dim U.
 (iv)  ODD: A^2 + A + I = 0 (so A is fixed-point-free of order 3, i.e. multiplication by omega on F_4^5,
       5 odd) and B != I.
 (v)   W0: rank r = 450, i.e. W(A, B) = 0.
Prints one summary with the counts per criterion; exit status 1 on any failure.
usage: python3 -I proofs/check_certificates.py [CERTFILE]
"""
import json, os, sys

n = 10
PAIRS = [(i, j) for j in range(n) for i in range(j)]
NP = len(PAIRS)                     # 45
NU = NP * n                         # 450
PIDX = {p: q for q, p in enumerate(PAIRS)}
ID = [1 << i for i in range(n)]


# ---------------------------------------------------------------- GF(2) linear algebra on ints
def rank(vecs):
    """Rank of a list of GF(2) vectors (ints)."""
    piv = {}                        # leading bit -> vector
    r = 0
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v; r += 1
                break
    return r


def in_span(v, vecs):
    return rank(list(vecs) + [v]) == rank(vecs)


def mat_vec(M, v):                  # M given by columns
    r = 0
    for i in range(n):
        if (v >> i) & 1:
            r ^= M[i]
    return r


def mat_mul(M, K):                  # columns of M K
    return [mat_vec(M, c) for c in K]


def mat_add(*Ms):                  # columns of the sum
    out = []
    for cs in zip(*Ms):
        r = 0
        for c in cs:
            r ^= c
        out.append(r)
    return out


def kernel_dim(M):                  # dim ker M = n - rank of the columns
    return n - rank(M)


def parity(x):
    return bin(x).count("1") & 1


# ---------------------------------------------------------------- the system defining W(A, B)
def system(A, B):
    """Rows (450-bit ints) of the linear system beta(A e_i, A e_j) + B beta(e_i, e_j) = 0, i < j."""
    eqs = []
    for (i, j) in PAIRS:
        x, y = A[i], A[j]
        # beta(x, y) = sum_{k<l} (x_k y_l + x_l y_k) beta(e_k, e_l)
        coef = [(((x >> k) & (y >> l)) ^ ((x >> l) & (y >> k))) & 1 for (k, l) in PAIRS]
        q0 = PIDX[(i, j)]
        for s in range(n):          # output bit s
            row = 0
            for q in range(NP):
                if coef[q]:
                    row ^= 1 << (q * n + s)
            for t in range(n):      # (B beta_ij)_s = sum_t B[t]_s * bit t of beta_ij
                if (B[t] >> s) & 1:
                    row ^= 1 << (q0 * n + t)
            eqs.append(row)
    return eqs


def pack(vals):                     # 45 pair values -> 450-bit vector
    v = 0
    for q, x in enumerate(vals):
        assert 0 <= x < (1 << n)
        v |= x << (q * n)
    return v


def form(vals):                     # beta as an n x n table of values beta(e_i, e_j)
    T = [[0] * n for _ in range(n)]
    for q, (i, j) in enumerate(PAIRS):
        T[i][j] = T[j][i] = vals[q]
    return T


def bil(T, x, y):
    r = 0
    for i in range(n):
        if (x >> i) & 1:
            for j in range(n):
                if (y >> j) & 1:
                    r ^= T[i][j]
    return r


# ---------------------------------------------------------------- checks
def check(c, types):
    tid = c["type"]; A = c["A"]; B = c["B"]; p = c["p"]; crit = c["criterion"]
    errs = []
    ok = lambda cond, msg: cond or errs.append(msg)

    ok(len(A) == n and len(B) == n and all(0 <= x < 1 << n for x in A + B), "A, B not 10 x 10")
    ok(tid in types and types[tid] == (p, A, B), "A, B, p differ from types10.tsv")
    ok(rank(A) == n and rank(B) == n, "A or B singular")
    Ap, Bp = ID, ID
    for _ in range(p):
        Ap, Bp = mat_mul(A, Ap), mat_mul(B, Bp)
    ok(Ap == ID and Bp == ID and (A != ID or B != ID), "(A, B) not of order p")

    # (i) full solution space
    eqs = system(A, B)
    r = rank(eqs)
    Wv = c["W"]
    ok(all(len(w) == NP for w in Wv), "W vectors must have 45 pair values")
    Wp = [pack(w) for w in Wv]
    ok(all(parity(e & w) == 0 for e in eqs for w in Wp), "a W vector does not solve the system")
    ok(rank(Wp) == len(Wp) == NU - r, "W basis: %d vectors of rank %d, expected 450 - %d = %d"
       % (len(Wp), rank(Wp), r, NU - r))
    Ts = [form(w) for w in Wv]

    if crit == "W0":
        ok(r == NU and not Wv, "W0: rank %d != 450" % r)
    elif crit in ("C", "WONG"):
        a, D = c["a"], c["D"]
        ok(0 < a < 1 << n and D and all(0 < d < 1 << n for d in D), "a or D out of range")
        ok(rank(D) == len(D), "D not linearly independent")
        if crit == "C" or c.get("invariant"):
            ok(all(in_span(mat_vec(A, d), D) for d in D), "D not A-invariant")
        E = rank([bil(T, a, d) for T in Ts for d in D])
        need = len(D) - (1 if in_span(a, D) else 0)
        ok(E < need, "compression fails: dim E = %d, dim D - [a in D] = %d" % (E, need))
        ok(E == c["dim_E"], "stored dim_E differs")
    elif crit == "R":
        D = c["D"]; k = len(D)
        ok(rank(D) == k and k >= 3, "D not independent or dim D < 3")
        ok(all(in_span(mat_vec(A, d), D) for d in D), "D not A-invariant")
        E = rank([bil(T, D[i], D[j]) for T in Ts for j in range(k) for i in range(j)])
        ok(E <= k - 1, "restriction fails: dim E_W(D) = %d > k - 1 = %d" % (E, k - 1))
        ok(E == c["dim_E"], "stored dim_E differs")
    elif crit == "Q":
        ok(p == 3, "Q needs p = 3")
        U = kernel_dim(mat_add(mat_mul(A, A), A, ID))
        FB = kernel_dim(mat_add(B, ID))
        ok(2 * FB < U, "order three fails: dim U = %d, dim Fix B = %d" % (U, FB))
        ok((U, FB) == (c["dim_U"], c["dim_FixB"]), "stored dims differ")
    elif crit == "ODD":
        ok(mat_add(mat_mul(A, A), A, ID) == [0] * n, "A^2 + A + I != 0")
        ok(kernel_dim(mat_add(A, ID)) == 0, "A has fixed points")
        ok(B != ID, "B = I")
        ok(n % 4 == 2, "n / 2 not odd")
    else:
        errs.append("unknown criterion %r" % crit)
    return errs


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "certificates10.jsonl")
    types = {}
    for l in open(os.path.join(here, "..", "types", "types10.tsv")):
        f = l.rstrip("\n").split("\t")
        types[int(f[0])] = (int(f[1]), list(map(int, f[5].split(","))), list(map(int, f[6].split(","))))
    order = ("W0", "C", "R", "Q", "WONG", "ODD")
    cnt = dict.fromkeys(order, 0); bad = 0; seen = set()
    for line in open(path):
        c = json.loads(line)
        errs = check(c, types)
        if c["type"] in seen:
            errs.append("duplicate type")
        seen.add(c["type"])
        if errs:
            bad += 1
            print("FAIL type %s (%s): %s" % (c["type"], c["criterion"], "; ".join(errs)), flush=True)
        else:
            cnt[c["criterion"]] = cnt.get(c["criterion"], 0) + 1
    total = sum(cnt.values())
    print("certificates checked: %d types (%d of the %d types), %d failed" % (total + bad, len(seen), len(types), bad))
    print("first criterion: " + ", ".join("%s %d" % (k, cnt[k]) for k in order)
          + "   (Table III: W=0 31, C 237, R 4, Q 4, Wong 2, odd 2; total 280)")
    expected = {"W0": 31, "C": 237, "R": 4, "Q": 4, "WONG": 2, "ODD": 2}
    if bad or cnt != expected:
        print("CERTIFICATE CHECK FAILED" if bad else "counts differ from Table III")
        sys.exit(1)
    print("ALL %d CERTIFICATES VALID" % total)


if __name__ == "__main__":
    main()
