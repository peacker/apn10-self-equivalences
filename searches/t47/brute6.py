#!/usr/bin/env python3
"""Brute-force validation of the type-47 solver in dimension six (paper, Section "Type 47", last paragraph).

Setting (the n = 6 analogue, as in gold6.py): A = B = x -> x^4 on F_64 (order 3), in the basis
b1, A b1, b2, A b2, 1, g^21 (bits 0..5). U = ker(A^2 + A + 1) = F_4^2 (bits 0..3, F_4-coordinate k = bits 2k, 2k+1
as 1, omega; f_0 = bit 0, f_1 = bit 2 is an F_4-basis of U), Z = Fix(A) = F_4 (bits 4, 5: z_0, z_1).
The U x U part (kappa, h) of beta is fixed to that of x^3 (file g6_3.txt = `gold6.py 3`, pair index -> value).
The remaining blocks of an equivariant beta are
  tau: U x Z -> U' (F_4-linear in u):  E[j][k] = tau(f_k, z_j), 4 entries of 4 bits = 16 bits
  Q:   Z x Z -> Z' (alternating):      q = Q(z_0, z_1), 2 bits
(U x Z -> Z' and Z x Z -> U' vanish by equivariance). This script runs over all 2^18 values of (E, q), builds the
LUT of each F and tests APN. Then it applies the normalization of the solver:
  - E_0 = (tau(f_0, z_0), tau(f_0, z_1)) is the reduced echelon basis of X = tau(f_0, Z) (c0 < c1, different
    leading bits, c1 zero at the leading bit of c0);
  - X is the representative of its orbit under the F_4-linear output maps P' with P' kappa = kappa (the one with
    the smallest bitset of its 4 elements),
and compares the normalized APN completions with the solver output g6_found.txt (stage2g / stage2i at n = 6, see
verify.py). Expected (2026-11-02 run, reconstructed from the session log): 144 APN completions, 24 with E_0 in
echelon form, 6 with X a representative, equal to the 6 LUTs of g6_found.txt; all 6 have the ortho-derivative
spectrum of the Gold function x^3.
usage: python3 -I searches/t47/brute6.py [g6_3.txt [g6_found.txt]]   (a few minutes; exit status 1 on mismatch)"""
import os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
uu_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "g6_3.txt")
found_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "g6_found.txt")

# F_4 = {0, 1, omega, omega^2} as 2-bit integers (1, omega), omega^2 = omega + 1
M4 = [[0] * 4 for _ in range(4)]
for a in range(4):
    for b in range(4):
        r = 0
        for i in range(2):
            if b >> i & 1:
                r ^= a << i
        if r & 4:
            r ^= 7
        M4[a][b] = r


def coord(u, k):                 # F_4-coordinate k of u in U
    return (u >> (2 * k)) & 3


def sm(c, v):                    # scalar c in F_4 times v in U
    return M4[c][coord(v, 0)] | (M4[c][coord(v, 1)] << 2)


b = [[0] * 4 for _ in range(4)]  # U x U part of beta (pair index p = j(j-1)/2 + i)
for line in open(uu_file):
    p, v = map(int, line.split())
    j = 1
    while j * (j + 1) // 2 <= p:
        j += 1
    i = p - j * (j - 1) // 2
    b[i][j] = b[j][i] = v


def bUU(x, y):
    r = 0
    for i in range(4):
        if x >> i & 1:
            for j in range(4):
                if y >> j & 1:
                    r ^= b[i][j]
    return r


def lut_of(E, q):
    """LUT of F with F(v) = sum_{i<j} v_i v_j beta(e_i, e_j)."""
    beta = [[0] * 6 for _ in range(6)]
    for i in range(4):
        for j in range(4):
            beta[i][j] = b[i][j]
    for i in range(4):           # tau(e_i, z_j), e_i in U
        for j in range(2):
            t = sm(coord(1 << i, 0), E[j][0]) ^ sm(coord(1 << i, 1), E[j][1])
            beta[i][4 + j] = beta[4 + j][i] = t
    beta[4][5] = beta[5][4] = q << 4
    F = [0] * 64
    for v in range(64):
        r = 0
        for j in range(6):
            if v >> j & 1:
                for i in range(j):
                    if v >> i & 1:
                        r ^= beta[i][j]
        F[v] = r
    return F


def is_apn(F):
    return all(len({F[x] ^ F[x ^ a] for x in range(64)}) == 32 for a in range(1, 64))


allapn = 0
echelon = {}
for code in range(1 << 18):
    E = [[(code >> (4 * (2 * j + k))) & 15 for k in range(2)] for j in range(2)]   # E[j][k] = tau(f_k, z_j)
    q = code >> 16
    F = lut_of(E, q)
    if not is_apn(F):
        continue
    allapn += 1
    c0, c1 = E[0][0], E[1][0]
    if c0 and c1 and c0 < c1 and c0.bit_length() != c1.bit_length() and not (c1 >> (c0.bit_length() - 1) & 1):
        echelon[",".join(map(str, F))] = (c0, c1)

# output maps: F_4-linear P' on U' (P'(f_0) = p0, P'(f_1) = p1), invertible, fixing kappa = U'-part of bUU
GO = []
for p0 in range(16):
    for p1 in range(16):
        P = lambda v, p0=p0, p1=p1: sm(coord(v, 0), p0) ^ sm(coord(v, 1), p1)
        if len({P(v) for v in range(16)}) < 16:
            continue
        if all(P(b[i][j] & 15) == b[i][j] & 15 for i in range(4) for j in range(4)):
            GO.append(P)


def sset(c0, c1):                # bitset of X = span{c0, c1}
    return sum(1 << v for v in {0, c0, c1, c0 ^ c1})


rep = {s for s, (c0, c1) in echelon.items() if all(sset(P(c0), P(c1)) >= sset(c0, c1) for P in GO)}
found = {l.strip() for l in open(found_file) if l.strip()}


# Gold check: ortho-derivative differential spectrum equal to that of x^3 on F_64 (an EA-invariant)
def ortho_spec(L):
    pi = [0] * 64
    for a in range(1, 64):
        img = {L[x] ^ L[x ^ a] ^ L[a] ^ L[0] for x in range(64)}
        pi[a] = next(c for c in range(1, 64) if all(bin(c & v).count("1") % 2 == 0 for v in img))
    C = Counter()
    for a in range(1, 64):
        cnt = Counter(pi[x] ^ pi[x ^ a] for x in range(64))
        for y in range(64):
            C[cnt.get(y, 0)] += 1
    return sorted(C.items())


def gmul(x, y):
    r = 0
    while y:
        if y & 1:
            r ^= x
        y >>= 1; x <<= 1
        if x & 64:
            x ^= 0b1011011       # x^6 + x^4 + x^3 + x + 1
    return r


gold = [gmul(x, gmul(x, x)) for x in range(64)]
sg = ortho_spec(gold)
all_gold = all(ortho_spec(list(map(int, s.split(",")))) == sg for s in rep)

print("APN completions of the 2^18:", allapn)
print("with E_0 in reduced echelon form:", len(echelon))
print("|output stabilizer of kappa|:", len(GO), "; with X an orbit representative:", len(rep))
print("equal to the solver output (g6_found.txt, %d LUTs):" % len(found), rep == found)
print("all of them Gold (ortho-derivative spectrum of x^3):", all_gold)
ok = rep == found and all_gold
print("OK" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
