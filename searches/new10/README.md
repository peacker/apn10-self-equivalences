# Candidate new quadratic APN classes at n = 10 (type 73, found 2026-10-29)

Type 73: order 5, A = B = identity on 2 coordinates + two 4x4 blocks with minimal polynomial x^4+x^3+x^2+x+1
(conjugate to the Frobenius x -> x^4, so these are the quadratic functions with coefficients in F_4)
(server run, sym10b/out4/t_73.*.out, full setting, about 5% explored on 2026-10-29).

| fingerprint | file line | M | permutation test | ortho-derivative differential spectrum |
|---|---|---|---|---|
| cc9135ebccd0fd0d | t73.luts line 1 (key K0fc8b78adad833ff) | 3 | none (M < 5) | 0:635846 2:317181 4:79085 6:13540 8:1725 10:145 12:30 |
| dcffd811adb7bb59 | t73.luts line 2 (key Kdc011a5b5ff63d6d) | 4 | none (M < 5) | 0:633636 2:320446 4:79035 6:12670 8:1580 10:180 12:5 |

- Checked: differential uniformity 2, algebraic degree 2 (independent code), invariants differ from all 25 known
  classes (apndb/n10/quadratic_classes.tsv) and from each other, so they are 2 classes outside those 25 (proved:
  the invariant is a CCZ invariant).
- Not CCZ-equivalent to a permutation (counting criterion, M < 5).
- Open: literature after our survey (paper E4/E11) may contain them. Do not claim or publish before that check.
- Files: t73.luts (LUTs), t73.raw / t73_all.raw (raw engine output), spectra.py.

## Univariate form (univariate.py, modulus x^10 + x^3 + 1, F_4 = {0, 1, 236, 237})
After conjugating A to the Frobenius x -> x^4, both functions have all coefficients in F_4 (output of
univariate.py; the linear terms can be dropped up to EA-equivalence). This is the "(10, 2)" case of
Ichanska-Berg-Kaleyski-Yu (ePrint 2024/1838), whose exhaustive search they report as infeasible (partial runs
only, no new class reported). It is Beierle-Leander's class 140 (code github.com/cbe90/quadratic_apn): their
randomized search there found nothing in about 16 CPU-days.

## Literature check (agent, 2026-10-29)
No source up to October 2026 lists an n = 10 quadratic APN class beyond our 25; no printed spectrum matches.
Not readable (403): Helleseth-Ichanska-Kaleyski ePrint 2026/2401 (switching families that "generalize previously
known sporadic APN instances in dimension 10": check by hand), Kaleyski-Sunde 2026/940, Yu et al. 2024/1778.

## Helleseth-Ichanska-Kaleyski ePrint 2026/2401 (read 2026-10-30)
Their n = 10 functions F_l, F_1, F_2, F_3 (Sec. IV, hik26.py) are triplicate #1, #6, #4, #3 of our 25; their G_1,
G_2 are BHK (F13-type, in our 25). Our spectra equal their Table IV. Neither candidate matches. Both added to
apndb (commit 97efc43); to be included as original results in our paper (PAPER_NOTES.md).
Note: these two LUTs satisfy F(Ax) = A F(x) only up to a linear term (engine output convention); univariate.py now
prints the quadratic part only. Quadratic parts: 32 and 35 terms, all coefficients in F_4.

## Type 73 complete (2026-11-02)
The re-split run (192 pieces, server run sym10b/out5) finished on 2026-11-02 07:55 UTC (T73_192_DONE; every piece
has its "done:" line). 4,536 leaves, 24 distinct engine keys, 24 classes by the ortho-derivative invariant
(`sym10/t73_final/classes.txt`): 20 of the 25 published classes and our 4 new ones (cc9135eb, dcffd811, 0a0ea0bf,
fc409cf6). So the quadratic APN functions of F_2^10 with coefficients in F_4 form exactly 24 EA-classes up to
our invariant (it separates all 29 known classes). None is CCZ-equivalent to a permutation (the new ones have
M < 5; the known ones are known not to be).
EA check (2026-11-02, t73_final/ea_check.txt): the 515 printed leaves form exactly 24 EA-classes, each equivalent
to its apndb representative (certificates re-checked). The other 4,021 leaves were only fingerprinted (same
ortho-derivative spectra as these 24). See PAPER_NOTES.md for the exact claim.
