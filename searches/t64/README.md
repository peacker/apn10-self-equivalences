# Type 64: exploratory searches (not exhaustive)

Type 64 (A = multiplication by omega on F_4^5, B = I; the "canonical triplicates") is **open**. Nothing in this
directory is a classification. The searches are exploratory; they give evidence for the conjecture of the paper
(every APN function of type 64 is EA-equivalent to one of the 16 known classes), not a proof.

## What is checked exactly

- `iso64.py` (in `verify.py`): exactly 16 of the 29 known classes have a self-equivalence of type 64 (an isometry
  of order 3 without fixed points), none of the 4 new ones.
- `herm.py`: the Hermitian criterion (APN iff every point of PG(4, 4) is killed by exactly one nonzero element of
  V) on x^3, x^9 and x^3 + Tr(x^9).
- `rankdist64.py` (in `verify.py`): rank distributions of V and of C = V^perp for the 16 known classes.

## What the searches establish, and what they do not

All searches below keep, for speed, only the **first** APN output per *fingerprint* (the differential spectrum of
the ortho-derivative, a cheap invariant; `sw2.c`, `tri10/tri10.c`, `poly64.c`) and discard the others. The kept
outputs are then named by the full ortho-derivative invariant (`classify_sw2.py`, `classify_poly.py`,
`tri10/tri10_job.py`) compared with the 29 known classes; **no EA-equivalence certificate is computed**.
The fingerprint is not a complete invariant. So the searches establish only:

> every APN function output by these searches has the fingerprint of one of the known classes of type 64
> (and the first output per fingerprint has the full invariant of that class).

A new class whose fingerprint equals that of a known class would not be detected. The searches are not
independent of each other in any formal sense, and none of them covers all of type 64 (the generic search tree has
a Knuth estimate of about 4.5e57 nodes).

| Search | Files | Scope | Outputs |
|---|---|---|---|
| Sparse polynomials, exhaustive within the listed shapes | `tri10/tri10.c`, `tri10_job.py`, `tri10_plan.sh` | x^b + sum c_k x^(e_k): b = 3 with c_k in F_4 and 5 <= t <= 10, b = 3 with c_k in F_32 and 1 <= t <= 4, b in {9, 33} with c_k in F_4 and 5 <= t <= 9 | `tri10/results/` (classes per shard, job counts) |
| Codimension-2 Hermitian switching | `sw2.c`, `run_sw2.sh`, `seeds64.py`/`seeds64.txt`, `classify_sw2.py` | one space V per known class (16 seeds); for each of the 174,251 8-dim subspaces V_8 of V, all 10-dim V' != V containing V_8 that satisfy the Hermitian criterion. **sw2 excludes V itself** (a completion of every V_8). This contains every codimension-1 switching. | `sw2out/` (`keys.tsv`: first space per fingerprint and its name) |
| Random sparse polynomials | `poly64.c`, `classify_poly.py` | random polynomials with 2 to 14 of the 25 type-64 monomials, coefficients in F_4, F_32 or F_1024 | `polyout/` |
| Codimension-1 switching BFS (superseded) | `switch64.py`, `bfs64.py` | class-level BFS from the 16 known classes; **capped**: at most `max_new = 8` new spaces per hyperplane, `reps = 2` representatives expanded per class, and the solutions modulo the hyperplane are enumerated from m = 1 upward and cut at the cap | `bfs64.out`, `sw64.out` |

The codimension-1 BFS is superseded by `sw2`, which computes, for each seed, all codimension-1 (and codimension-2)
switchings exactly, without a cap. Since EA-equivalence need not respect the F_4-structure, a class may contain
several inequivalent spaces V; `sw2` switched from one space per class only.

`t64.c` is a point-by-point exhaustive solver (each V containing H_1 = diag(0, 1, 1, 1, 1) exactly once); it was
used for Knuth estimates (`leaves_knuth.txt`) and cannot finish on type 64. `mc4.py` is a small Monte-Carlo
experiment on 4 x 4 Hermitian matrices (rank of random subspaces), not a search of type 64.
