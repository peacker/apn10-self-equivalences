#!/usr/bin/env python3
"""Affine search space for the self-equivalence type A = B = J_m on F_2^n (n = 2m; Fix(A) = Im N_A of dim m), with
the intertwiner reductions:
 (a) beta on Fix(A) x Fix(A) = the bilinear form of the Gold function x^(2^k+1) on F_{2^m} (every quadratic APN
     (m, m)-function is EA-equivalent to a Gold one for the m used here, and C(A), C(B) act on the fixed spaces as
     the full GL(m));
 (b) for x in Fix(A) \\ 0 and w_x with N_A w_x = x: beta(x, w_x) lies in Fix(B) and must avoid the hyperplane
     I_x = beta(x, Fix(A)) (otherwise T_x = beta(x, .) has a kernel larger than {0, x}): one affine equation per x;
 (c) normalisation: the unipotent symmetries w -> w + Q w (Q: complement -> Fix(A)) and y -> y + U y (U: complement
     of Fix(B) -> Fix(B)) translate the Fix(B)-components of beta(Fix(A), complement) by a known subspace L; the
     coordinates at the pivots of L are set to 0 (one representative per translation class: exact).
Input: the reordered-basis W-file of sym_setup.py (first m basis vectors span Fix(A)). Output: a W-file for sym10
(echelon basis of the homogeneous part, no normalisation data) and an offset file "pair value" (a particular
solution); run sym10 with -nonorm -offset OFFSET and no -symdepth / -cp.
usage: intertwine.py WFILE TYPES_FILE TYPE_ID n k MODULUS OUT_W OUT_OFFSET"""
import sys
wfile, tfile, tid, n, k, mod, out_w, out_off = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6], 0), sys.argv[7], sys.argv[8]
m = n // 2
N = 1 << n
pairs = [(i, j) for j in range(n) for i in range(j)]
pidx = {p: q for q, p in enumerate(pairs)}
P = len(pairs)

# --- W basis (beta as an integer: bit p * n + r) ---
lines = open(wfile).read().split("\n")
hd = lines[0].split(); assert int(hd[0]) == n
d = int(hd[1])
W = []
for t in range(d):
    f = lines[1 + t].split()
    vals = list(map(int, f[2:2 + P]))
    v = 0
    for q in range(P):
        v |= vals[q] << (q * n)
    W.append(v)
def pairval(beta, i, j):
    if i == j: return 0
    if i > j: i, j = j, i
    return (beta >> (pidx[(i, j)] * n)) & (N - 1)

# --- A' (reordered) and B ---
row = [l.rstrip("\n").split("\t") for l in open(tfile) if l.split("\t")[0] == tid][0]
B = list(map(int, row[6].split(",")))
sys.path.insert(0, __file__.rsplit("/", 1)[0] + "/../types")
import subprocess
Ap = list(map(int, subprocess.run([sys.executable, __file__.rsplit("/", 1)[0] + "/../types/reorder.py", tid, tfile, str(n)],
                                  capture_output=True, text=True).stdout.strip().split(",")))
def app(M, v):
    r = 0
    for i in range(n):
        if (v >> i) & 1: r ^= M[i]
    return r
for i in range(m): assert app(Ap, 1 << i) == 1 << i, "first m basis vectors must be fixed by A'"
NA = [Ap[i] ^ (1 << i) for i in range(n)]
NB = [B[i] ^ (1 << i) for i in range(n)]
def rank(vs):
    rows = list(vs); r = 0
    for b in range(n):
        pv = next((i for i in range(r, len(rows)) if (rows[i] >> b) & 1), None)
        if pv is None: continue
        rows[r], rows[pv] = rows[pv], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> b) & 1: rows[i] ^= rows[r]
        r += 1
    return r
fixB = []
for v in range(1, N):
    if app(B, v) == v and rank(fixB + [v]) > len(fixB): fixB.append(v)
assert len(fixB) == m
compB = []
for i in range(n):
    if rank(fixB + compB + [1 << i]) > len(fixB) + len(compB): compB.append(1 << i)
coordB = {}                                       # vector -> (fix coords, comp coords)
for a in range(1 << m):
    for c in range(1 << m):
        v = 0
        for t in range(m):
            if (a >> t) & 1: v ^= fixB[t]
            if (c >> t) & 1: v ^= compB[t]
        coordB[v] = (a, c)
def fixvec(a):
    v = 0
    for t in range(m):
        if (a >> t) & 1: v ^= fixB[t]
    return v

# --- Gold form on F_{2^m} with polynomial basis ---
F = 1 << m
def fmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a & F: a ^= mod
    return r
def fpow(a, e):
    r = 1
    while e:
        if e & 1: r = fmul(r, a)
        a = fmul(a, a); e >>= 1
    return r
qk = 1 << k
def G(x, y):                                       # x, y in F_{2^m} (coordinates of Fix(A)), value in Fix(B)
    return fixvec(fmul(fpow(x, qk), y) ^ fmul(x, fpow(y, qk)))

# --- linear system over t in F_2^d: rows (coefficient bitmask over t, rhs) ---
rows = []
def eq_bit(pq, r, rhs):                             # bit r of beta on pair index pq equals rhs
    c = 0
    for t in range(d):
        if (W[t] >> (pq * n + r)) & 1: c |= 1 << t
    rows.append((c, rhs))
for j in range(m):                                  # (a) Gold on Fix x Fix
    for i in range(j):
        g = G(1 << i, 1 << j)
        for r in range(n): eq_bit(pidx[(i, j)], r, (g >> r) & 1)
def solve(rows):
    piv = {}
    for c, rhs in rows:
        for b, (pc, pr) in piv.items():
            if (c >> b) & 1: c ^= pc; rhs ^= pr
        if not c:
            if rhs: return None
            continue
        b = c.bit_length() - 1
        for b2 in list(piv):
            pc, pr = piv[b2]
            if (pc >> b) & 1: piv[b2] = (pc ^ c, pr ^ rhs)
        piv[b] = (c, rhs)
    t0 = 0
    for b, (c, rhs) in piv.items():
        if rhs: t0 |= 1 << b
    free = [t for t in range(d) if t not in piv]
    null = []
    for f_ in free:
        v = 1 << f_
        for b, (c, rhs) in piv.items():
            if (c >> f_) & 1: v |= 1 << b
        null.append(v)
    return t0, null
def beta_of(t):
    v = 0
    for i in range(d):
        if (t >> i) & 1: v ^= W[i]
    return v
sol = solve(rows); assert sol, "Gold prefix inconsistent"
beta_ref = beta_of(sol[0])

# (b) hyperplane conditions for every x in Fix(A) \ 0
nb = 0
for x in range(1, F):
    # w_x in span(e_m..e_{n-1}) with N_A w_x = x (x in coordinates of e_0..e_{m-1})
    wx = next(v << m for v in range(1, 1 << (n - m)) if app(NA, v << m) == x)
    Ix = [G(x, 1 << i) for i in range(m)]
    rIx = rank(Ix)
    assert rIx == m - 1
    ell = next(l for l in range(1, N) if all(bin(l & u).count("1") % 2 == 0 for u in Ix)
               and any(bin(l & fixvec(a)).count("1") % 2 for a in range(1, F)))
    c = 0
    for t in range(d):
        val = 0
        for i in range(m):
            if (x >> i) & 1:
                for j in range(m, n):
                    if (wx >> j) & 1: val ^= pairval(W[t], i, j)
        if bin(ell & val).count("1") % 2: c |= 1 << t
    rows.append((c, 1)); nb += 1

# (c) normalisation: translations of the Fix(B)-parts of beta(e_i, e_j), i < m <= j
shifts = []
for a in range(m):                                  # Q_{a,b}: e_{m+b} -> e_a
    for b in range(n - m):
        v = 0
        for i in range(m):
            val = G(1 << i, 1 << a) if i != a else 0
            v |= val << (pidx[(i, m + b)] * n)
        shifts.append(v)
for a in range(m):                                  # U_{a,b}: compB[b] -> fixB[a]
    for b in range(m):
        v = 0
        for i in range(m):
            for j in range(m, n):
                fc, cc = coordB[pairval(beta_ref, i, j)]
                if (cc >> b) & 1: v |= fixB[a] << (pidx[(i, j)] * n)
        shifts.append(v)
ech = []                                            # reduced echelon of L (pivot = lowest bit)
for v in shifts:
    for e in ech:
        if v & (e & -e): v ^= e
    if v:
        lo = v & -v
        ech = [e ^ v if e & lo else e for e in ech]
        ech.append(v)
for e in ech:
    pos = (e & -e).bit_length() - 1
    eq_bit(pos // n, pos % n, 0)
sol = solve(rows); assert sol, "constraints inconsistent"
t0, null = sol
beta0 = beta_of(t0)
H = [beta_of(v) for v in null]
hech = []                                           # echelon in search order (pivot = lowest bit)
for v in H:
    for e in hech:
        if v & (e & -e): v ^= e
    if v:
        lo = v & -v
        hech = [e ^ v if e & lo else e for e in hech]
        hech.append(v)
hech.sort(key=lambda e: (e & -e).bit_length())
with open(out_w, "w") as f:
    f.write(f"{n} {len(hech)}\n")
    for e in hech:
        pv = (e & -e).bit_length() - 1
        f.write(f"{pv // n} {pv % n} " + " ".join(str((e >> (q * n)) & (N - 1)) for q in range(P)) + "\n")
    f.write("0\n0\n0\n")
with open(out_off, "w") as f:
    for q in range(P):
        f.write(f"{q} {(beta0 >> (q * n)) & (N - 1)}\n")
print(f"type {tid} (n = {n}, Gold k = {k}): W dim {d}; equations: Gold prefix {m * (m - 1) // 2 * n}, "
      f"hyperplane {nb}, normalisation {len(ech)}; remaining free bits {len(hech)}")
