#!/usr/bin/env python3
"""Proof-based exclusions of the n = 10 self-equivalence types, computed from the space W of equivariant forms.

For a type (A, B), W is the space of alternating bilinear maps beta with beta(Ax, Ay) = B beta(x, y) (sym_setup.py).
F is APN iff for every a != 0 the map beta(a, .) has kernel {0, a}. The criteria below only use W, so each
exclusion is a proof for every function of the type.

 (C) compression: for a != 0 and a subspace D, every beta in W maps a x D into S(a, D) = span{beta_w(a, d)}
     (w in a basis of W, d in a basis of D). APN needs dim D - [a in D] <= dim S(a, D). D runs over the sums of
     the invariant subspaces ker g(A) (g | minimal polynomial). This contains the intertwiner lemma (a in Fix A)
     and the dual lemma (a in an isotypic part, D = Fix A).
 (R) restriction + Chevalley-Warning: for an invariant D with beta(D, D) in E = span{beta_w(d, d')}: F restricted
     to D is an APN (k, dim E)-function, k = dim D. dim E = k - 1 is impossible for k >= 3: the Pfaffian of
     the alternating forms e^* beta|D (e^* in E^*, a (k-1)-dim space) has degree k/2 < k - 1 variables, so
     Chevalley-Warning gives a nonzero e^* with a degenerate form, i.e. a kernel larger than {0, a}.
     (For odd k every alternating form is degenerate, so already dim E >= k - 1 is too weak: see (C).)
 (Q) order 3: for x in U = ker(A^2 + A + 1), q(x) = beta(x, Ax) lies in Fix(B) and q(x) != 0 for x != 0; q is a
     system of dim Fix(B) quadratic forms in dim U variables, so Chevalley-Warning needs 2 dim Fix(B) >= dim U.
 (P) probe: max rank of beta_t(a, .) over random t in W; a value < n - 1 for some a flags a type that is
     probably empty (not a proof).
usage: proofs10.py WDIR [TYPES...]   (WDIR/w_ID.txt from sym_setup.py; prints one line per type)"""
import os, random, sys
import numpy as np

n = 10; N = 1 << n
here = os.path.dirname(os.path.abspath(__file__))
rows = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(os.path.join(here, "..", "types", "types10.tsv"))}


def app(M, v):
    r = 0
    for i in range(n):
        if (v >> i) & 1:
            r ^= M[i]
    return r


def span_basis(vs):
    basis = []
    for v in vs:
        v = int(v)
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return basis


def poly_cols(g, M):            # columns of g(M), g as an integer bit polynomial
    cols = []
    for i in range(n):
        p = 1 << i; acc = 0; k = 0
        while (g >> k):
            if (g >> k) & 1:
                acc ^= p
            p = app(M, p); k += 1
        cols.append(acc)
    return cols


def kernel(cols):               # basis of {v : sum v_i cols[i] = 0}
    aug = [(cols[i], 1 << i) for i in range(n)]
    basis = []                  # Gaussian elimination on (image, preimage)
    piv = []
    ker = []
    for img, pre in aug:
        for (bi, bp) in piv:
            if img ^ bi < img:
                img ^= bi; pre ^= bp
        if img:
            piv.append((img, pre))
            piv.sort(reverse=True)
        else:
            ker.append(pre)
    return span_basis(ker)


def pmod(a, b):
    db = b.bit_length()
    while a and a.bit_length() >= db:
        a ^= b << (a.bit_length() - db)
    return a


def irreducible_factors(f):     # factors (with multiplicity) of a squarefree-or-not polynomial over F2, small degree
    out = []
    g = 2
    while f.bit_length() > 1:
        if g.bit_length() > f.bit_length():
            out.append(f); break
        if pmod(f, g) == 0:
            # divide
            q = 0; a = f; db = g.bit_length()
            while a and a.bit_length() >= db:
                s = a.bit_length() - db; q |= 1 << s; a ^= g << s
            out.append(g); f = q
        else:
            g += 1
    return out


def invariant_family(A):
    """Sums of the kernels ker g(A)^k over the irreducible factors g of the minimal polynomial, k = 1, 2."""
    # the order of A is 1, 2 or an odd prime p: the factors of x^ord - 1 cover the minimal polynomial
    order = 1; P = list(A)
    while P != [1 << i for i in range(n)] and order < 2000:
        P = [app(A, v) for v in P]; order += 1
    gens = set()
    for g in sorted(set(irreducible_factors((1 << order) | 1))):
        for k in (1, 2):
            gk = 1
            for _ in range(k):          # gk = g^k
                r = 0; a, b = gk, g
                while b:
                    if b & 1: r ^= a
                    a <<= 1; b >>= 1
                gk = r
            K = kernel(poly_cols(gk, A))
            if K:
                gens.add(tuple(sorted(K)))
    gens = [list(g) for g in gens]
    fam = set()
    m = len(gens)
    for mask in range(1, 1 << m):
        vs = []
        for i in range(m):
            if (mask >> i) & 1:
                vs += gens[i]
        fam.add(tuple(sorted(span_basis(vs))))
    return [list(D) for D in fam]


def load_w(path):
    L = open(path).read().split()
    nn, d = int(L[0]), int(L[1]); assert nn == n
    P = n * (n - 1) // 2; pos = 2
    T = np.zeros((d, n, n), dtype=np.uint16)
    for w in range(d):
        pos += 2
        vals = L[pos:pos + P]; pos += P
        q = 0
        for j in range(n):
            for i in range(j):
                v = int(vals[q]); q += 1
                T[w, i, j] = v; T[w, j, i] = v
    return T


def apply_vec(T, a):             # (d, n): beta_w(a, e_j)
    M = np.zeros(T.shape[:1] + (n,), dtype=np.uint16)
    for i in range(n):
        if (a >> i) & 1:
            M ^= T[:, i, :]
    return M


def images(M, D):                # all beta_w(a, d) for d in the basis D
    out = []
    for v in D:
        x = np.zeros(M.shape[0], dtype=np.uint16)
        for j in range(n):
            if (v >> j) & 1:
                x ^= M[:, j]
        out.append(x)
    return np.unique(np.concatenate(out)) if out else np.zeros(0, dtype=np.uint16)


def in_span(a, D):
    return len(span_basis(D + [a])) == len(D)


def rank_rows(M):                 # rank of the n x n matrix whose columns are given as ints
    return len(span_basis(M))


def analyse(tid, wpath, probes=24):
    r = rows[tid]
    A = list(map(int, r[5].split(","))); B = list(map(int, r[6].split(",")))
    T = load_w(wpath)
    d = T.shape[0]
    found = []
    fam = invariant_family(A)
    if d == 0:
        return ["W=0"], None
    # (C) compression
    best = None
    for a in range(1, N):
        M = apply_vec(T, a)
        for D in fam:
            S = span_basis(images(M, D))
            need = len(D) - (1 if in_span(a, D) else 0)
            if len(S) < need:
                best = (a, D, len(S), need); break
        if best:
            break
    if best:
        a, D, s, need = best
        found.append("C(a=%d,dimD=%d,dimS=%d<%d)" % (a, len(D), s, need))
    # (R) restriction + Chevalley-Warning
    for D in fam:
        k = len(D)
        if k < 3:
            continue
        vals = []
        for x in D:
            Mx = apply_vec(T, x)
            vals.append(images(Mx, D))
        E = span_basis(np.unique(np.concatenate(vals)))
        if len(E) == k - 1:
            found.append("R(dimD=%d,dimE=%d)" % (k, len(E))); break
    # (Q) order 3
    if r[1] == "3":
        U = kernel(poly_cols(0b111, A)); FB = kernel(poly_cols(0b11, B))
        if U and 2 * len(FB) < len(U):
            found.append("Q(dimU=%d,dimFixB=%d)" % (len(U), len(FB)))
    # (P) probe of the commutative rank
    rng = random.Random(int(tid))
    ts = [np.array([rng.getrandbits(1) for _ in range(d)], dtype=bool) for _ in range(probes)]
    Bm = [np.bitwise_xor.reduce(T[t], axis=0) if t.any() else np.zeros((n, n), dtype=np.uint16) for t in ts]
    worst = n
    for a in range(1, N):
        mx = 0
        for Mt in Bm:
            cols = [0] * n
            for i in range(n):
                if (a >> i) & 1:
                    for j in range(n):
                        cols[j] ^= int(Mt[i, j])
            mx = max(mx, rank_rows(cols))
            if mx == n - 1:
                break
        worst = min(worst, mx)
        if worst < n - 1:
            break
    return found, worst


if __name__ == "__main__":
    wdir = sys.argv[1]
    ids = sys.argv[2:] or list(rows)
    for tid in ids:
        f, worst = analyse(tid, "%s/w_%s.txt" % (wdir, tid))
        r = rows[tid]
        print("%s\tp=%s\t%s\t%s\t%s\tprobe_maxrank_min=%s\t%s" % (tid, r[1], r[2], r[3], r[4], worst, " ".join(f) or "-"),
              flush=True)
