# Type 24: notes on the dedicated solver (`c/full24g.c`)

Type 24 is the set of quadratic functions of F_{2^10} with all coefficients in F_32. This note states what the
solver relies on; the proofs are in the paper (Section "Type 24", Lemmas "Normal form", "Planes", "regular").
F_32 uses the polynomial x^5 + x^2 + 1.

## Normal form

Write x = (s, t) with s, t in F_2^5, with A(s, t) = (s, s + t) and B = A on the output. Every quadratic APN
function with this self-equivalence is EA-equivalent to

    F(s, t) = ( G(s) + nu(s, s) + mu(s, t),  R(s) + nu(t, s) + G(t) ),

with G the Gold function x^3 or x^5 of F_32 (env `GE` = 3 or 5), mu the polar form of G, nu bilinear (125 bits)
and R quadratic with polar form rho (alternating, 50 bits). In the code: `x = p | q << M`, `value = top | bot << M`,
`S` is nu, `Dp` is rho (`form24.py` checks the normal form).

## APN as one condition per plane

For independent u, v let Sigma(u, v) be the set of nu(b, v) + nu(c, u) + mu(b, c) over all b, c with
mu(u, v) + nu(u, v) + nu(v, u) + mu(u, c) + mu(b, v) = 0, and for a plane P let Forb(P) be the union of
Sigma(u, v) over the three bases {u, v} of P. F is APN iff (i) nu(b, b) is not in Im mu(b, .) for every b != 0
(the "diag conditions") and (ii) rho(u, v) is not in Forb(P) for every plane P = <u, v> (155 planes).
Forb(P) depends only on the columns nu(., w), w in P.

## Translation space T (dim 49)

The maps (s, t) -> (s, t + Y s) on the input and (u, v) -> (u, v + Y' u) on the output (Y, Y' in M_5(F_2)) commute
with A, keep G and mu, change rho, and replace nu by nu + theta, theta(x, y) = mu(x, Y y) + Y' mu(x, y). These theta
span a space T of dimension 49 (`setup_T` computes a reduced echelon basis). nu is taken as the unique element of
nu + T that vanishes on the pivot coordinates.
Invariance: under a translation, (b, c) -> (b + Y u, c + Y v) is a bijection of the solutions in the definition of
Sigma(u, v), and Sigma(u, v) is shifted by a vector that depends only on the plane <u, v>. So each part of Forb(P),
and Forb(P), is translated by one vector: sizes, intersections and the affine dimension do not change.

## Automorphisms of G

The maps L(x) = lambda x^(2^k) (lambda in F_32^*, 0 <= k < 5) satisfy G o L = N o G with N linear. Such a map
sends F to F' = diag(N^{-1}, N^{-1}) o F o diag(L, L) (inverse N on both output halves, so that mu is kept:
N^{-1} mu(L x, L y) = mu(x, y)). It commutes with A and replaces nu by nu'(x, y) = N^{-1} nu(L x, L y). Hence
Forb_F'(P) = N^{-1} Forb_F(L P), with the three parts permuted, and key_F'(P) = key_F(L P).

Regular action (paper, Lemma "regular"). The 155 maps act regularly on the 155 planes of F_32. Proof: it suffices
that only the identity fixes a plane P. The group has order 155, prime to 3! = 6, so an element g fixing P (it
permutes the 3 nonzero points of P) fixes each of them. If k = 0, lambda a = a with a != 0 gives lambda = 1. If
k != 0, the nonzero fixed points satisfy x^(2^k - 1) = lambda^(-1); since gcd(2^k - 1, 31) = 1 this has exactly one
solution, so g has at most one nonzero fixed point, a contradiction. `planes_check.py` recomputes that the 155 maps
are automorphisms of both x^3 and x^5 and send <e_0, e_1> to 155 distinct planes (all of them).

## The plane key

For a plane P the key is (|Forb(P)|, sorted sizes of its three parts), then, to break ties, the affine dimension,
the sizes of the pairwise and triple intersections of the parts and a hash of their multiplicity histograms
(`ckey`, `pkey`). By the two sections above it is invariant under translations of nu and transported by the
automorphisms. The production run uses `MINKEY=1` (see `run24.sh`): the solver keeps only prefixes in which no plane
inside the current span has a strictly smaller key than <e_0, e_1> (ties are kept).

## Coverage

Let F be APN of type 24. By the normal form it is EA-equivalent to a function as above. Choose a plane P_0 of
smallest key. By the regular action there is exactly one automorphism L of G with L <e_0, e_1> = P_0. Then
key_F'(<e_0, e_1>) = key_F(P_0) <= key_F(L P) = key_F'(P) for every plane P. Then reduce nu modulo T;
this does not change any key. The result has <e_0, e_1> of smallest key among all planes, so it passes every key
test (which rejects only strictly smaller keys), and it satisfies (i) and (ii), so the exhaustive search below the
key tests reaches it as a leaf. Hence every class of type 24 has at least one leaf. The run is split into 4,096 jobs
(GE = 3, 5; PART = r/2048: the column-1 options with index r mod 2048); `check_run.py` checks that exactly these jobs
finished with exit code 0, and lists the 3 leaves.
