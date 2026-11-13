#!/usr/bin/env python3
"""Independent check of sw2.c (codimension-2 Hermitian switching): for a sample of the 8-dim subspaces U of a seed
space V, count all 10-dim V' containing U that satisfy the Hermitian criterion, by a different method, and compare
with the number of spaces V' != V that sw2 reports for the same subspaces U.

Nothing is imported from sw2.c or switch64.py; only the encoding of Hermitian matrices (25-bit codes, as in sw2.c and
herm.py) and of F_4^5 (2 bits per coordinate, 0, 1, w, w^2 = 0, 1, 2, 3) is shared.
Criterion: V' is APN iff, for every point q of PG(4, 4), the map H -> Hq from V' to q^perp has rank 9.
Every completion V' != V found here is also checked to be APN with the full difference table.
Method: choose a complement K of U (17-dim), reduce Kx q modulo Uq for all 2^17 vectors x of K at once; the plane
<a, b> of K gives a good V' iff for every q the reductions of aq, bq span a space of dimension 9 - dim Uq.

usage: check_sw2.py SW2_BINARY SEED_LINE N STEP START   (seed line of seeds64.txt, N subspaces U with index = START
mod STEP, in the order of sw2; exit 1 on a mismatch)"""
import itertools, os, re, subprocess, sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
MUL = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
CONJ = [0, 1, 3, 2]


def herm_matrix(code):
    """5 x 5 matrix over F_4 of a 25-bit code: bits 0-4 diagonal, then for i < j the entries (1, 1) and (w, w^2)."""
    M = [[0] * 5 for _ in range(5)]
    pos = [(b, b, 1, 1) for b in range(5)]
    for i, j in itertools.combinations(range(5), 2):
        pos += [(i, j, 1, 1), (i, j, 2, 3)]
    for b, (i, j, u, v) in enumerate(pos):
        if code >> b & 1:
            if i == j:
                M[i][i] ^= u
            else:
                M[i][j] ^= u
                M[j][i] ^= v
    for i in range(5):
        for j in range(5):
            assert M[j][i] == CONJ[M[i][j]]
    return M


def apply(M, q):
    """M q for q in F_4^5 (2 bits per coordinate) -> 10-bit integer."""
    c = [(q >> 2 * j) & 3 for j in range(5)]
    r = 0
    for i in range(5):
        s = 0
        for j in range(5):
            s ^= MUL[M[i][j]][c[j]]
        r |= s << 2 * i
    return r


def echelon(vecs):
    """echelon basis (dict leading bit -> vector) of the span of integer bit vectors"""
    basis = {}
    for v in vecs:
        for p in sorted(basis, reverse=True):
            if v >> p & 1:
                v ^= basis[p]
        if v:
            p = v.bit_length() - 1
            for p2 in list(basis):
                if basis[p2] >> p & 1:
                    basis[p2] ^= v
            basis[p] = v
    return basis


def reduce_all(arr, basis):
    """reduce every entry of arr modulo the span of an echelon basis (fully reduced, so the result is canonical)"""
    arr = arr.copy()
    for p in sorted(basis, reverse=True):
        hit = (arr >> p) & 1 == 1
        arr[hit] ^= basis[p]
    return arr


def subspaces_in_sw2_order():
    """the 2-dim subspaces {f1, f2, f1 ^ f2} of the dual of V, in the order of sw2.c"""
    for f1 in range(1, 1024):
        for f2 in range(f1 + 1, 1024):
            if f1 ^ f2 >= f2:
                yield f1, f2


def count_good(V, f1, f2, points, imgs):
    par = lambda x: bin(x).count("1") & 1
    T = [t for t in range(1, 1024) if not par(f1 & t) and not par(f2 & t)]
    Ub = list(echelon(T).values())
    assert len(Ub) == 8
    U = []
    for t in Ub:
        h = 0
        for k in range(10):
            if t >> k & 1:
                h ^= V[k]
        U.append(h)
    EU = echelon(U)
    assert len(EU) == 8
    Kb = []                                           # complement of U: greedy over unit vectors, low bits first
    for b in range(25):
        if len(echelon(U + Kb + [1 << b])) > 8 + len(Kb):
            Kb.append(1 << b)
    assert len(Kb) == 17
    n = 1 << 17
    zero_d1 = []
    inC = np.ones(n, dtype=bool)
    inC[0] = False
    for q in range(len(points)):
        Uq = echelon([imgs[q](h) for h in U])
        d = 9 - len(Uq)
        arr = np.zeros(1, dtype=np.int64)
        for v in Kb:                                  # Kx q for all x, x indexed by its 17 coordinates
            arr = np.concatenate([arr, arr ^ imgs[q](v)])
        red = reduce_all(arr, Uq)
        assert len(echelon(set(red.tolist()) | set(Uq.values()))) == 9   # everything lies in q^perp (dim 9)
        if d == 2:
            inC &= red != 0
        elif d == 1:
            zero_d1.append(red == 0)
        elif d != 0:
            return 0, None
    Z = np.packbits(np.array(zero_d1), axis=0).T if zero_d1 else None    # shape (n, bytes)
    C = np.nonzero(inC)[0]
    good = 0
    planes = []
    for a in C:
        b = C[C > a]
        c = b ^ a
        b = b[(c > b) & inC[c]]
        if len(b) and Z is not None:
            b = b[~(Z[b] & Z[a]).any(axis=1)]
        good += len(b)
        planes += [(int(a), int(x)) for x in b]
    # the plane that gives V itself: both of its vectors lie in V
    def vec(x):
        h = 0
        for i in range(17):
            if x >> i & 1:
                h ^= Kb[i]
        return h
    vplanes = 0
    for a, b in planes:
        if all(len(echelon(V + [vec(x)])) == 10 for x in (a, b)):
            vplanes += 1
        else:                                         # every other completion: APN by the full difference table
            assert is_apn(U + [vec(a), vec(b)]), (f1, f2, a, b)
    return good, vplanes


def is_apn(W):
    """F(x) = (x^* H_k x)_k for the 10 Hermitian matrices of W; APN iff F(x + a) + F(x) takes 512 values for all a"""
    Ms = [herm_matrix(h) for h in W]
    F = np.zeros(1024, dtype=np.int64)
    for x in range(1024):
        v = 0
        for k, M in enumerate(Ms):
            y = apply(M, x)
            s = 0
            for i in range(5):
                s ^= MUL[CONJ[(x >> 2 * i) & 3]][(y >> 2 * i) & 3]
            assert s < 2
            v |= s << k
        F[x] = v
    xs = np.arange(1024)
    return all(len(np.unique(F[xs ^ a] ^ F)) == 512 for a in range(1, 1024))


def main():
    exe, line, N, step, start = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    V = list(map(int, open(os.path.join(here, "seeds64.txt")).read().splitlines()[line - 1].split("\t")[1].split()))
    points = [v for v in range(1, 1024) if [(v >> 2 * i) & 3 for i in range(5) if (v >> 2 * i) & 3][0] == 1]
    assert len(points) == 341
    cache = {}

    def img(q):
        def f(h):
            key = (q, h)
            if key not in cache:
                cache[key] = apply(herm_matrix(h), points[q])
            return cache[key]
        return f
    imgs = [img(q) for q in range(341)]
    total = 0
    picked = 0
    for idx, (f1, f2) in enumerate(subspaces_in_sw2_order()):
        if idx >= N * step:
            break
        if idx % step != start % step:
            continue
        g, vp = count_good(V, f1, f2, points, imgs)
        assert vp == 1, (idx, vp)                     # V itself is a good completion of every U, exactly once
        total += g - 1
        picked += 1
    r = subprocess.run([exe] + list(map(str, V)) + ["-umax", str(N), "-ustep", str(step), "-ustart", str(start)],
                       capture_output=True, text=True)
    m = re.search(r"subspaces (\d+), new spaces (\d+)", r.stderr)
    sw = int(m.group(2))
    ok = sw == total
    print(f"seed {line}: {picked} subspaces U; independent count of V' != V: {total}; sw2: {sw}; "
          f"{'agree' if ok else 'DIFFER'}", flush=True)
    sys.exit(0 if ok else 1)


main()
