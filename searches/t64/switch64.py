#!/usr/bin/env python3
"""Search for new classes in type 64 by Hermitian switching.

SUPERSEDED, NOT USED AS EVIDENCE IN THE PAPER. As run (bfs64.py), neighbours() stops after max_new solutions per
hyperplane U and enumerates the solution coset from small combinations upward, so the run was not exhaustive; only
reps representatives per class were expanded. sw2.c does every codimension-1 (and -2) switching exactly, without a
cap (independent check: check_sw2.py). This file is kept for lut_of and inv_hash, used by classify_*.py.

An APN function of type 64 is a 10-dim space V of Hermitian 5x5 matrices over F_4 such that every point of PG(4, 4)
is killed (H q = 0) by exactly one nonzero element of V. Codimension-1 switching: keep a hyperplane U of V and
replace the tenth basis element by H'. Points killed by U need H' q not in U q, a hyperplane of q^perp (9-dim):
one affine equation l_q(H' q) = 1. Points not killed by U are killed exactly once by the new coset automatically.
So all neighbours sharing U are the solutions of one affine system (exact only without the max_new cap).
Classes are named by the ortho-derivative invariant (classify.py's hash) and compared with apndb/n10.
usage: switch64.py [--rounds R] [--reps K] OUT_PREFIX   (BFS over classes; K representatives expanded per class)"""
import hashlib, os, random, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402

F4M = [[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]; CONJ = [0, 1, 3, 2]
def coord(v, i): return (v >> (2 * i)) & 3

# points of PG(4, 4): first nonzero coordinate 1
PTS = []
for v in range(1, 1024):
    lead = next(coord(v, i) for i in range(5) if coord(v, i))
    if lead == 1: PTS.append(v)
assert len(PTS) == 341

def basis_mat(b):
    M = [[0] * 5 for _ in range(5)]
    if b < 5: M[b][b] = 1; return M
    k = 5
    for i in range(5):
        for j in range(i + 1, 5):
            if b == k: M[i][j] = 1; M[j][i] = 1; return M
            if b == k + 1: M[i][j] = 2; M[j][i] = 3; return M
            k += 2
def matvec(M, v):
    r = 0
    for i in range(5):
        s = 0
        for j in range(5): s ^= F4M[M[i][j]][coord(v, j)]
        r |= s << (2 * i)
    return r
BM = [basis_mat(b) for b in range(25)]
IMG = [[matvec(BM[b], q) for b in range(25)] for q in PTS]      # IMG[q][b] = E_b q (10-bit)
def hq(H, qi):
    r = 0; b = 0
    while H:
        if H & 1: r ^= IMG[qi][b]
        H >>= 1; b += 1
    return r

def lut_of(V):
    """F(x) = (x^* H_k x)_k on F_4^5 (10-bit index as in PTS encoding)."""
    # quadratic form of basis matrix b at x
    q = np.zeros((25, 1024), dtype=np.uint8)
    xs = np.arange(1024)
    c = [(xs >> (2 * i)) & 3 for i in range(5)]
    f4 = np.array(F4M, dtype=np.uint8); cj = np.array(CONJ, dtype=np.uint8)
    tr = np.array([0, 0, 1, 1], dtype=np.uint8)
    for b in range(25):
        if b < 5:
            q[b] = f4[c[b], cj[c[b]]] & 1           # x_i conj(x_i) in F_2
            continue
        k = 5
        for i in range(5):
            for j in range(i + 1, 5):
                for h, bb in ((1, k), (2, k + 1)):
                    if bb == b:
                        q[b] = tr[f4[cj[c[i]], f4[h, c[j]]]]   # Tr(conj(x_i) h x_j)
                k += 2
    F = np.zeros(1024, dtype=np.int64)
    for kk, H in enumerate(V):
        bit = np.zeros(1024, dtype=np.uint8)
        for b in range(25):
            if (H >> b) & 1: bit ^= q[b]
        F |= bit.astype(np.int64) << kk
    return F

def inv_hash(F):
    return hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]

def span_basis(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b.append(v)
    return b

def solve_affine(rows, rhs, nbits=25):
    """All solutions of r . x = c (over F_2): returns (particular, basis of homogeneous) or None."""
    piv = {}                                   # leading bit -> (row, rhs)
    for r, c in zip(rows, rhs):
        for lb in sorted(piv, reverse=True):
            if (r >> lb) & 1: pr, pc = piv[lb]; r ^= pr; c ^= pc
        if r:
            lb = r.bit_length() - 1
            for k in list(piv):
                pr, pc = piv[k]
                if (pr >> lb) & 1: piv[k] = (pr ^ r, pc ^ c)
            piv[lb] = (r, c)
        elif c:
            return None
    x0 = 0
    for lb, (r, c) in piv.items():
        if c: x0 |= 1 << lb
    free = [b for b in range(nbits) if b not in piv]
    hom = []
    for f in free:
        v = 1 << f
        for lb, (r, c) in piv.items():
            if (r >> f) & 1: v |= 1 << lb
        hom.append(v)
    return x0, hom

def kill_vector(V, qi):
    """coordinates t (10-bit) of the unique nonzero element of V killing point qi (0 if none)."""
    imgs = [hq(H, qi) for H in V]
    for t in range(1, 1 << len(V)):
        s = 0
        for k in range(len(V)):
            if (t >> k) & 1: s ^= imgs[k]
        if s == 0: return t
    return 0

def neighbours(V, max_new=64):
    """All V' sharing a hyperplane with V (up to max_new per hyperplane), as lists of 10 matrices."""
    imgs = [[hq(H, qi) for H in V] for qi in range(341)]
    tq = []
    for qi in range(341):
        t = 0
        for cand in range(1, 1024):
            s = 0
            for k in range(10):
                if (cand >> k) & 1: s ^= imgs[qi][k]
            if s == 0: t = cand; break
        tq.append(t)
    out = []
    for f in range(1, 1024):                   # U = {t : popcount(f & t) even}
        # basis of U in V-coordinates
        Ub_t = [t for t in (1 << k for k in range(10))]
        ub = []
        piv = next(k for k in range(10) if (f >> k) & 1)
        for k in range(10):
            if k == piv: continue
            t = 1 << k
            if bin(f & t).count("1") & 1: t ^= 1 << piv
            ub.append(t)
        def elem(t):
            s = 0
            for k in range(10):
                if (t >> k) & 1: s ^= V[k]
            return s
        U = [elem(t) for t in ub]
        rows = []; rhs = []
        for qi in range(341):
            if bin(f & tq[qi]).count("1") & 1: continue          # H_q not in U: point not killed by U
            Uq = span_basis([hq(H, qi) for H in U])               # 8-dim
            # functional l on F_2^10 with l(Uq) = 0 and l(q^perp) != 0: kernel(l) = Uq + <something outside q^perp>
            # q^perp = {y : Tr(q^* y) = 0}; pick y0 outside q^perp; l = annihilator of Uq + <y0>
            qv = PTS[qi]
            def trp(y):
                s = 0
                for i in range(5): s ^= F4M[CONJ[coord(qv, i)]][coord(y, i)]
                return 1 if s in (2, 3) else 0
            y0 = next(y for y in range(1, 1024) if trp(y))
            Kb = span_basis(Uq + [y0])                             # 9-dim
            # annihilator: l with l . k = 0 for k in Kb (10-bit functionals)
            l = next(l for l in range(1, 1024) if all(bin(l & k).count("1") % 2 == 0 for k in Kb))
            # row: bit b of H' contributes l . IMG[q][b]
            r = 0
            for b in range(25):
                if bin(l & IMG[qi][b]).count("1") & 1: r |= 1 << b
            rows.append(r); rhs.append(1)
        sol = solve_affine(rows, rhs)
        if sol is None: continue
        x0, hom = sol
        # solutions modulo U: the homogeneous space contains U; extra directions give new neighbours
        Ubasis = span_basis(U)
        extra = []
        for h in hom:
            x = h
            for u in sorted(Ubasis, reverse=True): x = min(x, x ^ u)
            if x and x not in extra:
                tmp = span_basis(Ubasis + extra + [x])
                if len(tmp) > len(Ubasis) + len(extra): extra.append(x)
        if not extra: continue
        cnt = 0
        for m in range(1, 1 << len(extra)):
            Hn = x0
            for k in range(len(extra)):
                if (m >> k) & 1: Hn ^= extra[k]
            out.append(U + [Hn]); cnt += 1
            if cnt >= max_new: break
    return out

if __name__ == "__main__":
    args = sys.argv[1:]
    rounds = 3; reps = 2
    if "--rounds" in args: i = args.index("--rounds"); rounds = int(args[i + 1]); del args[i:i + 2]
    if "--reps" in args: i = args.index("--reps"); reps = int(args[i + 1]); del args[i:i + 2]
    out = args[0] if args else os.path.join(here, "sw64")
    known = {}
    for line in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")):
        t = line.rstrip("\n").split("\t")
        if t[0] != "name": known[t[1]] = t[0]
    sys.path.insert(0, here)
    from herm import code_of, gpow, gmul, n as _n          # seeds from GF(2^10)
    def trf(z):
        s = 0; x = z
        for _ in range(10): s ^= x; x = gmul(x, x)
        return s
    # seeds as V in our coordinates: rebuild from the F_4-basis used in herm.py
    import herm
    seeds = []
    for F in (lambda x: gpow(x, 3), lambda x: gpow(x, 9), lambda x: gpow(x, 3) ^ trf(gpow(x, 9))):
        V, _ = code_of(F)
        # herm.py uses the same 25-bit layout (diagonal bits, then (i, j) pairs with 2 bits: value bits b0, b1)
        seeds.append(V)
    classes = {}                                   # inv -> list of V
    for V in seeds:
        h = inv_hash(lut_of(V)); classes.setdefault(h, []).append(V)
    frontier = list(classes)
    log = open(out + ".log", "a")
    for rnd in range(rounds):
        newfront = []
        for h in frontier:
            for V in classes[h][:reps]:
                for W in neighbours(V):
                    F = lut_of(W)
                    hw = inv_hash(F)
                    if hw not in classes:
                        classes[hw] = [W]; newfront.append(hw)
                        name = known.get(hw, "NEW")
                        msg = "round %d class %s %s V=%s" % (rnd, hw, name, " ".join(map(str, W)))
                        print(msg, flush=True); log.write(msg + "\n"); log.flush()
                        if name == "NEW":
                            with open(out + ".new_luts", "a") as fo: fo.write(",".join(map(str, F.tolist())) + "\n")
                    elif len(classes[hw]) < reps:
                        classes[hw].append(W)
        frontier = newfront
        print("round %d: %d classes (%d new this round)" % (rnd, len(classes), len(newfront)), flush=True)
        if not frontier: break
