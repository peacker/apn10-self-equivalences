"""Explicit EA-equivalence of quadratic APN functions, with a checkable certificate.

If G = L1 o F o L2 + A with L1, L2 linear permutations and A affine, the ortho-derivatives satisfy
pi_G = M o pi_F o L2 with M = (L1^T)^-1 (with respect to the dot product of integers as bit
vectors). We enumerate the pairs (L2, M) with this property by a depth-first search over the
images of a basis under L2, pruned with point invariants of the ortho-derivatives (the number of
preimages of pi(a) and the multiset of row a of the DDT of pi), and accept the first pair for
which G + L1 o F o L2 is affine. The search is exhaustive, so `None` means that F and G are
EA-inequivalent; for quadratic APN functions EA- and CCZ-equivalence coincide (Yoshiara).
Every returned pair (L1, L2) is re-checked by `check_certificate`, which does not rely on the
search.
"""
import numpy as np

from .invariants import ortho_derivative


def apply(cols, v):
    """Images of the integer vectors v under the linear map with column images cols."""
    v = np.asarray(v, dtype=np.int64)
    out = np.zeros_like(v)
    for i, c in enumerate(cols):
        out ^= ((v >> i) & 1) * int(c)
    return out


def rank(cols):
    rows, r = [int(c) for c in cols], 0
    for bit in range(max(rows + [0]).bit_length() - 1, -1, -1):
        p = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def is_affine(H, n):
    """H(x + y) + H(x) + H(y) + H(0) = 0 for all x and all basis vectors y."""
    X = np.arange(1 << n)
    return all(np.all(H[X ^ (1 << i)] ^ H ^ H[1 << i] ^ H[0] == 0) for i in range(n))


def check_certificate(F, G, n, L1, L2):
    """True iff L1, L2 are invertible and G + L1 o F o L2 is affine."""
    F, G = np.asarray(F, dtype=np.int64), np.asarray(G, dtype=np.int64)
    if rank(L1) != n or rank(L2) != n:
        return False
    X = np.arange(1 << n)
    return is_affine(G ^ apply(L1, F[apply(L2, X)]), n)


def _point_invariants(pi, n):
    N = 1 << n
    X = np.arange(N)
    pre = np.bincount(pi, minlength=N)[pi]
    out = []
    for a in range(N):
        v, c = np.unique(np.bincount(pi[X ^ a] ^ pi, minlength=N), return_counts=True)
        out.append((int(pre[a]),) + tuple(zip(v.tolist(), c.tolist())))
    return out


def kernel_partners(F, n):
    """P[x] = set of nonzero y != x lying in a common kernel of a component bilinear form with x.

    If G = L1 o F o L2 + affine, the kernel of component v of G is L2^-1 of the kernel of
    component L1^T v of F, so L2 maps P_G[y] onto P_F[L2 y]. For a Gold function the nonzero
    kernels are the lines of a Desarguesian spread, and P[x] has 2 elements.
    """
    N = 1 << n
    X = np.arange(N)
    beta = np.stack([F[X ^ (1 << j)] ^ F ^ F[1 << j] ^ F[0] for j in range(n)], axis=1)   # beta(x, e_j)
    P = [set() for _ in range(N)]
    for v in range(1, N):
        t = beta & v
        for s in (1, 2, 4, 8, 16, 32):
            t = t ^ (t >> s)
        ker = np.nonzero(np.all((t & 1) == 0, axis=1))[0]
        ker = [int(x) for x in ker if x]
        if len(ker) > 1:
            for x in ker:
                P[x].update(y for y in ker if y != x)
    return P


class _PartialLinearMap:
    """A linear map over GF(2) known on a subspace, given by reduced (source, image) pairs."""

    def __init__(self, rows=None):
        self.rows = dict(rows or {})

    def reduce(self, s, t=0):
        for p in sorted(self.rows, reverse=True):
            if s >> p & 1:
                a, b = self.rows[p]
                s, t = s ^ a, t ^ b
        return s, t

    def add(self, s, t):
        """Impose s -> t; False if this contradicts linearity."""
        s, t = self.reduce(s, t)
        if s == 0:
            return t == 0
        self.rows[s.bit_length() - 1] = (s, t)
        return True


def _transpose_inverse(cols, n):
    """Columns of (M^T)^-1 for the matrix M with column images cols, or None if singular."""
    A = [[(cols[j] >> i) & 1 for j in range(n)] for i in range(n)]       # M[i][j]
    T = [[A[j][i] for j in range(n)] + [int(i == k) for k in range(n)] for i in range(n)]
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if T[i][c]), None)
        if p is None:
            return None
        T[r], T[p] = T[p], T[r]
        for i in range(n):
            if i != r and T[i][c]:
                T[i] = [x ^ y for x, y in zip(T[i], T[r])]
        r += 1
    return [sum(T[i][n + j] << i for i in range(n)) for j in range(n)]


_CACHE = {}


def ea_equivalence(F, G, n, node_limit=10**7, transitive=False, kernels=False):
    """(L1, L2) with G + L1 o F o L2 affine, or None if F and G are EA-inequivalent.

    F and G must be quadratic APN functions given by their tables. Raises RuntimeError if the
    node limit is reached.

    transitive=True: the caller asserts that the linear parts L2 of the EA-automorphisms of F
    (F o L2 = L1 o F + affine) act transitively on the nonzero inputs, as for any representative
    of a Gold class (x -> a x). If (L1, L2) is a solution and A2 is such an automorphism with
    linear output part A1, then (L1 A1^-1, A2 L2) is a solution too. So the image of the first
    basis vector can be fixed to one point, and the search stays exhaustive.

    kernels=True: also require that L2 maps the kernel partners of G (see kernel_partners) onto
    those of F, and take the next basis vector among the partners of points already mapped when
    possible (its image is then one of the partners of the image). This is an exact constraint.
    """
    F, G = np.asarray(F, dtype=np.int64), np.asarray(G, dtype=np.int64)
    key = (n, F.tobytes())
    if key not in _CACHE:                     # F is usually a fixed class representative: compute its data once
        if len(_CACHE) > 64:
            _CACHE.clear()
        pF0 = ortho_derivative(F, n)
        _CACHE[key] = {"p": pF0, "i": _point_invariants(pF0, n)}
    cF = _CACHE[key]
    pF, pG = cF["p"], ortho_derivative(G, n)
    iF, iG = cF["i"], _point_invariants(pG, n)
    if sorted(iF[1:]) != sorted(iG[1:]):
        return None
    if transitive and len(set(iF[1:])) != 1:
        raise ValueError("transitive=True, but the nonzero points of F have different invariants")
    by_inv = {}
    for a in range(1, 1 << n):
        by_inv.setdefault(iF[a], []).append(a)
    X = np.arange(1 << n)
    nodes = 0
    if kernels:
        if "P" not in cF:
            cF["P"] = kernel_partners(F, n)
        PF, PG = cF["P"], kernel_partners(G, n)
        if sorted(len(p) for p in PF) != sorted(len(p) for p in PG):
            return None

    def next_basis(L):
        """(b, candidate images) for the next basis vector of the source."""
        if kernels:
            for y in sorted(L):
                for z in sorted(PG[y]):
                    if z not in L:
                        return z, sorted(c for c in PF[L[y]] if iF[c] == iG[z])
        b = next(1 << i for i in range(n) if (1 << i) not in L)
        return b, by_inv.get(iG[b], [])

    def dfs(k, L, M):
        nonlocal nodes
        nodes += 1
        if nodes > node_limit:
            raise RuntimeError("node limit reached")
        if k == n:
            L2 = [L[1 << i] for i in range(n)]       # L is defined on the whole space
            red = [M.reduce(1 << i) for i in range(n)]
            if any(s for s, _ in red):
                raise RuntimeError("the image of the ortho-derivative does not span F_2^n")
            L1 = _transpose_inverse([t for _, t in red], n)
            if L1 is None or not is_affine(G ^ apply(L1, F[apply(L2, X)]), n):
                return None
            return L1, L2
        b, cands = next_basis(L)
        images = set(L.values())
        if transitive and k == 0:
            cands = cands[:1]
        for c in cands:
            if c in images:
                continue
            M2, new, ok = _PartialLinearMap(M.rows), {}, True
            for x, lx in L.items():
                y, ly = x ^ b, lx ^ c
                if ly == 0 or iF[ly] != iG[y] or not M2.add(int(pF[ly]), int(pG[y])):
                    ok = False
                    break
                new[y] = ly
            if ok and kernels:
                full = {**L, **new}
                ok = all(full[z] in PF[ly] for y, ly in new.items() for z in PG[y] if z in full)
            if ok:
                r = dfs(k + 1, {**L, **new}, M2)
                if r:
                    return r
        return None

    return dfs(0, {0: 0}, _PartialLinearMap())
