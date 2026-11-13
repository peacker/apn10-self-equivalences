# Quadratic APN functions of F_2^10 with a prime-order self-equivalence: code and data

Code, generated data and stored results of the paper

> E. Bellini, *Quadratic APN functions of F_2^10 with a prime-order self-equivalence* (title provisional),
> in preparation for IEEE Transactions on Information Theory.

The paper classifies the quadratic APN functions F of F_2^10 that have a linear self-equivalence (A, B) of prime
order, i.e. beta(Ax, Ay) = B beta(x, y) for the polar form beta of F. Up to conjugation and powers there are 299
such types (A, B). Each type is settled by a proof criterion, by an exhaustive search, or by both. Type 64 stays
open. The repository contains everything needed to re-check each claim. `verify.py` re-runs the cheap parts and
checks the stored outputs of the expensive searches in about 6 minutes (see "Three ways to check"). The production searches took about
36,000 core-hours and are not re-run by `verify.py`.

## Requirements

- Python 3.10 or later with NumPy (`pip install -r requirements.txt`)
- a C compiler (`cc`; the C core of `lib/apnsearch` is compiled on first import into `~/.cache/apnsearch`,
  override with `APNSEARCH_BUILD=<dir>`)

## Environment

- Field: F_2^10 = F_{2^10} with the polynomial x^10 + x^3 + 1 (`0x409`), polynomial basis; vectors are integers.
- F_4 inside F_{2^10}: {0, 1, 236, 237}, with omega = 236.
- Type 24 (coefficients in F_32): F_32 with the polynomial x^5 + x^2 + 1.
- Searches: three Linux servers (Docker containers, Ubuntu 22.04) with 2 x AMD EPYC 7763 (256 threads, 1 TB RAM),
  C compiled with gcc `-O3 -march=native`, run at `nice 19` (see `engine/polite.sh`). The servers were shared; our
  container could use at most 200 threads per server (`cpu.max`). Every job is single-threaded; a core-hour is one
  hour of wall-clock time of one job on one hardware thread (`secs` in the `.done` files), not a dedicated core. The runs did not log the
  gcc version. Queried on 2026-11-13, all three servers report gcc `(Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0`,
  Ubuntu 22.04.5 LTS, Python 3.10.12 and NumPy 1.24.3; a binary compiled on one of them during the campaign
  records the same gcc in its `.comment` section.
- Local checks (tested versions): macOS 26.3 (Apple M5 Pro; Python is an x86_64 build running under Rosetta 2,
  so `cc` called from it builds x86_64), Apple clang 16.0.0 (clang-1600.0.26.6), Python 3.13.3, NumPy 2.2.4.
- Every sanitized subprocess environment (`verify.py`, `searches/t47/*.py`) sets `TMPDIR` to a temporary directory
  of the run, so the compiler's temporary files stay there (needed in sandboxes that forbid the system default).

## Three ways to check

The checks fall into three modes. They differ in what they establish and in cost.

**(a) Regression and stored-log consistency** (`python3 verify.py`, about 6 minutes on one core;
`--skip-ea` skips section 6). The cheap parts are re-run (types, W spaces, proof criteria, master table, the
engine on small types and on stored validation runs, the dedicated solvers at n = 6 and on planted inputs), and the
stored outputs of the expensive searches are checked for completeness and consistency: exact sets of job ids, each
record present once, exit code 0, no node-limit marker, agreement between reruns, and equality between the leaves
a search emitted and the leaves its EA records cover. This mode trusts the stored logs. It does **not** re-run the
production searches, and it re-checks EA-equivalence only on a sample (section 6: the 3 type-24 leaves, one leaf
per type-73 class, a few type-231 leaves). The EA results of all other leaves are read from the stored records.

**(b) Independent checks of stored certificates and representatives** (seconds to a minute each, all also run by
`verify.py`, mostly in its sections 2, 5 and 7). Each recomputes a claim from stored data with its own code:

```sh
python3 -I proofs/check_certificates.py   # the 280 exclusion certificates (Table III); imports nothing from the repo, 2 s
python3 -I types/recount_table2.py        # recount of Table II (299 types)
python3 -I searches/t24/kim6_f8.py        # n = 6: a x^5 + a x^33 + x^40 (a = g^54 in F_8) is EA-equivalent to Kim's function
python3 -I searches/t24/planes_check.py   # the 155 automorphisms of x^3 and x^5 act regularly on the planes of F_32
python3 -I searches/t47/brute6.py         # n = 6: type-47 solver against a brute force over all 2^18 completions, 15 s
python3 -I searches/t47/regen_stage1.py   # the 806 type-47 stage-1 classes regenerated and compared with cls/, 40 s
python3 -I searches/new10/recompute_M.py  # M of the 4 new classes, exhaustively (3, 4, 4, 4), 1 min
python3 -I tools/check_manifest.py        # SHA-256 of every tracked file (MANIFEST.sha256); or: shasum -a 256 -c MANIFEST.sha256
```

`verify.py` also recomputes the full ortho-derivative invariants (not only their hashes) of the 29 representatives
of `data/quadratic_classes.tsv`, checks the stated spectra of N1-N4, and checks that the 29 invariants are pairwise
distinct.

**(c) Full regeneration** of the searches and of the EA-equivalences (not part of `verify.py`; Linux servers;
see "Building and running the searches" and "Cost of the searches"): the generic `sym10` runs over the searched
types (about 3,000 core-hours in total), type 24 (`run24.sh`, about 7,500), type 47 (`run47n.sh`/`run47i.sh`,
about 13,300), type 73 (about 12,600, plus 5,800 for the rerun printing every leaf), the fingerprint-gap rerun with
an EA test of every leaf (`searches/fpgap/run_all.sh`: search 9, EA checks about 1,300 core-hours), and the EA test
of every type-73 leaf (`searches/t73/rerun_all/ea_all.py`, about 30 minutes).

`verify.py` prints the sanitized environment of its subprocesses, the Python, NumPy and compiler versions, one line
per check and `ALL CHECKS PASSED` at the end (exit status 1 otherwise). It unpacks the `.tar.gz` archives of stored
outputs next to themselves on first use (git ignores the unpacked copies).
Last run (2026-11-13): all 53 checks passed in 370 s (macOS 26.3, Apple M5 Pro, one core; archives already
unpacked).

### The certificate file

`proofs/certificates10.jsonl` has one line per type excluded by proof: the type, the criterion (the first that
applies, in the order of Table III: `W0`, `C` compression on invariant subspaces, `R` restriction, `Q` order three,
`WONG` compression with Wong sequences, `ODD` odd dimension over F_4), the matrices A and B (column bitmasks,
column i = A e_i), a basis of W(A, B) (per basis vector the 45 values beta(e_i, e_j), pairs ordered (0,1), (0,2),
(1,2), (0,3), ...) and the witness: a and a basis of D for C and WONG, a basis of D for R, dim U and dim Fix(B)
for Q. It is written by `proofs/make_certificates.py` (from `W/`, which `types/make_w.sh` writes; it uses the
search code of `proofs10.py` and `wong10.py`). `proofs/check_certificates.py` imports nothing from the
repository. For every certificate it rebuilds the linear system beta(Ax, Ay) = B beta(x, y) (450 unknowns) from
A and B, computes its rank r, checks that the stored basis has 450 - r independent solutions (so it spans all of
W(A, B)), and then checks the inequality of the criterion. Expected output:

```
certificates checked: 280 types (280 of the 299 types), 0 failed
first criterion: W0 31, C 237, R 4, Q 4, WONG 2, ODD 2   (Table III: W=0 31, C 237, R 4, Q 4, Wong 2, odd 2; total 280)
ALL 280 CERTIFICATES VALID
```

## Claims and how to re-check them

| Claim (paper) | Files | How to re-check | Cost |
|---|---|---|---|
| There are 299 types (A, B) of prime order at n = 10 (60 at n = 6, 107 at n = 8) | `types/enumerate.py`, `types/types{6,8,10}.tsv` | `python3 types/enumerate.py 10` (in `verify.py`) | seconds |
| The space W of polar forms of each type | `types/sym_setup.py`, `types/make_w.sh`; `engine/W3/` (the 35 spaces used by the searches) | `types/make_w.sh` writes `W/w_ID.txt` for all 299 types; deterministic, equal to `engine/W3` (in `verify.py`) | 1 min |
| Proof criteria C (compression), R (restriction + Chevalley-Warning), Q (order 3), W=0 prove 276 types | `proofs/proofs10.py` -> `proofs/proofs10.txt` | `python3 proofs/proofs10.py W` (in `verify.py`) | 10 s |
| Exclusion certificates of the 280 types excluded by proof (Table III: 31 W=0, 237 C, 4 R, 4 Q, 2 Wong, 2 odd) | `proofs/make_certificates.py` -> `proofs/certificates10.jsonl`; checker `proofs/check_certificates.py` | `python3 -I proofs/check_certificates.py` (in `verify.py`) | 2 s |
| Wong criterion proves types 130, 169 | `proofs/wong10.py` -> `proofs/wong10.txt` | `python3 proofs/wong10.py W 130 169` (in `verify.py`) | 1 s |
| Odd-dimension lemma ("lemma56") excludes types 55-59, in particular 55, 56 | `proofs/lemma56/check56.py` (numerical check of the lemma) | `python3 proofs/lemma56/check56.py 56` (in `verify.py`) | 1 s |
| Dual intertwiner lemma D excludes 95 types | `proofs/lemma_dual.py` -> `proofs/lemma_dual.txt` | in `verify.py` | 1 s |
| The bound of Baudrin et al. (DCC 2025, Cor. 3) excludes the 65 types with dim Fix(A) + dim Fix(B) >= 12 (among them type 41) | `proofs/bbf.py` -> `proofs/bbf_excluded.txt` | in `verify.py` | 1 s |
| Master table: 242 types by proof and search, 38 by proof only, 18 by search only, type 64 open | `proofs/master_table.py` (inputs: `types/types10.tsv`, `proofs10.txt`, `lemma_dual.txt`, `bbf_excluded.txt`, `search_status.txt`, `server_done.txt`) -> `proofs/master_table.tsv` | in `verify.py` | 1 s |
| The generic search engine is exhaustive and reproducible | `engine/sym10.c`, `engine/sym10f73.c`, `engine/validation/` | type 231 re-run (91,636 nodes, 961 leaves, as on the server); 12 forced type-73 runs (in `verify.py`) | 20 s |
| Search-only types with APN functions (129, 170, 171, 175, 176, 180, 231): every leaf printed again and EA-tested | `searches/fpgap/` | `python3 searches/fpgap/check_fp.py` (in `verify.py`); full re-run: `run_all.sh` | rerun complete (2026-11-08): 67 of 67 jobs agree, all 4,851,099 leaves EA-equivalent to a known class; search 9 core-hours, EA checks 1,332 core-hours; check 2 s |
| **Type 24** (coefficients in F_32): exactly 3 classes, Gold x^3, x^9 and x^3 + Tr(x^9) | `searches/t24/` (normal form `form24.py`, plane lemma, solver `c/full24g.c`, `run24.sh`, outputs `results/t24_out.tar.gz`, EA check `results/ea_check24.txt`) | `python3 searches/t24/check_run.py`; EA certificates recomputed in `verify.py`; full re-run: 4,096 jobs of `run24.sh` | 7,500 core-hours; check 1 min |
| **Type 47**: empty | `searches/t47/` (stage 1 `stage1.py`, `stage1b.py` -> 806 classes `cls/`, regenerated by `regen_stage1.py`; rank-h lemma; solvers `stage2g.c`, `stage2i.c`; outputs `results/t47_out.tar.gz`) | `python3 searches/t47/check_coverage.py` (in `verify.py`): 703 classes dead, 8,056 slices of the other 103 all searched, no leaf | 13,261 core-hours; check 1 min |
| **Type 73** (coefficients in F_4): exactly 24 classes, 20 published and 4 new | `searches/t73/final/` (first run, 192 pieces), `searches/t73/rerun_all/` (rerun printing every leaf, 3,072 shards; `ea_all.txt`: every one of the 4,536 leaves EA-tested) | `check73.py` (exact piece and shard id sets; rerun = first run, piece by piece; leaf set = leaf set of `ea_all.txt`, in `verify.py`); one EA certificate per class recomputed; full EA re-check: `python3 searches/t73/rerun_all/ea_all.py` (about 21 min; rerun from an empty `ea_all.txt` on 2026-11-13: all 4,536 leaves EA-equivalent to their class, certificates checked, 24 classes) | 12,600 + 5,800 (rerun) core-hours |
| The 4 new classes are APN, quadratic, have M = 3, 4, 4, 4, and differ from all 25 published classes | `data/quadratic_classes.tsv` (copy of apndb `n10/quadratic_classes.tsv`), `searches/new10/` (univariate forms over F_4; shortest forms in `short_forms.txt`, checked by `check_short_forms.py`; `recompute_M.py`) | in `verify.py` (APN test, degree, ortho-derivative spectra, full invariants of all 29 pairwise distinct, M recomputed) | 1 min |
| Type 41 is empty | BBF (`proofs/bbf.py`), criteria C and D; the partial search is in `searches/t41/` (not needed) | in `verify.py` | - |
| **Type 64** (canonical triplicates) is open; 16 of the 29 known classes lie in it; exploratory searches (not exhaustive, first output per fingerprint only, no EA certificates) output only functions with the fingerprints of known classes; see `searches/t64/README.md` | `searches/t64/` (`iso64.py`, `herm.py`, `t64.c`, `switch64.py`, `bfs64.py`, `sw2.c`, `poly64.c`, `tri10/`) | `python3 searches/t64/iso64.py` (in `verify.py`); `python3 searches/t64/herm.py` (Hermitian criterion on x^3, x^9, x^3 + Tr(x^9)); `python3 searches/t64/rankdist64.py` (rank distributions of V and C for the 16 classes, in `verify.py`); `check_sw2.py` (independent count of the codimension-2 switchings on sampled subspaces, in `verify.py`) | 3 min |

## Cost of the searches

The heavy searches took about 36,000 core-hours (column `core_hours` of `proofs/master_table.tsv`: 36,378
core-hours, from the `.done` files of the finished runs). Per type:

| Type | Search | Core-hours |
|---|---|---|
| 47 | dedicated solver (stage2g, stage2i), 8,056 slices | 13,261 |
| 73 | `sym10`, 192 pieces (plus 381 for the first 32-piece run; the rerun printing every leaf took about 5,800 more) | 12,600 |
| 24 | dedicated solver `full24g`, 4,096 jobs | 7,500 |
| 98 | `sym10` | 1,385 |
| 93 | `sym10` | 953 |
| 100 | `sym10` | 335 |
| the other 254 searched types | `sym10` | 344 in total |

For type 47 we count the summed run times of the stored finished jobs. Earlier stage2g work on slices that was
stopped and redone with stage2i (about 3,100 more core-hours) is not counted. Also not counted above: the fingerprint-gap rerun (`searches/fpgap`), the
partial type-41 runs, the type-64 searches and the Knuth estimates.

## Layout

| Path | Content |
|---|---|
| `types/` | `enumerate.py` (types of prime-order (A, B), Beierle-Leander style), `recount_table2.py` (independent recount of Table II), `types{6,8,10}.tsv`, `sym_setup.py` (the space W and the centraliser bases, input of `sym10`), `make_w.sh`, `reorder.py`, `classify.py` (name LUTs by the ortho-derivative invariant) |
| `proofs/` | `proofs10.py`/`.txt`, `wong10.py`/`.txt`, `make_certificates.py`, `certificates10.jsonl`, `check_certificates.py`, `lemma_dual.py`/`.txt`, `lemma56/check56.py`, `bbf.py`/`bbf_excluded.txt`, `master_table.py`/`.tsv`, `search_status.txt` (per type: which server runs are complete), `server_done.txt` (first line of every server `.done` file; machine names replaced by `srv1`, `srv2`) |
| `engine/` | `sym10.c` (generic symmetric DFS over W; usage in the header), `sym10f73.c` (same search, faster checkpoints, `-all` prints every leaf), run scripts (`run_type.sh`, `run_fast.sh`, `run_cfg.sh`, `run_strong.sh`, `run_both.sh`, `polite.sh`), Knuth estimates (`costs10.tsv`, `knuth_both.tsv`, `knuth_sym.tsv`, ...), `W3/` (W spaces of the hard types), batch files, `intertwine.py`, `gen_force.py`, `force_*.txt`, `validation/` (runs on known types, n = 6, 8 and forced type-73 prefixes) |
| `searches/t24/` | type 24: `NOTES.md` (normal form, translations, plane key, regular action, coverage), `planes_check.py`, `kim6_f8.py` (n = 6: Kim's class has a representative with coefficients in F_8), `form24.py`, `topbent.py`, `sparse_sub.py`, `c/full24.c` (n = 6 validation), `c/full24g.c` (production), `run24.sh`, `check_run.py`, `results/` |
| `searches/t47/` | type 47: `stage1.py`, `stage1b.py`, `cls/` (806 classes), `s2/` (rank-h filter), `stage2.c`, `stage2g.c`, `stage2i.c`, `gold6.py` and `brute6.py` (n = 6 validation), `regen_stage1.py`, job lists, run scripts, `check_coverage.py`, `results/` |
| `searches/t73/` | type 73: `final/` (first complete run: `classes.txt`, `ea_check.py`/`.txt`, LUTs, outputs `out5.tar.gz`), `rerun_all/` (rerun with every leaf: `run73a.sh`, `check73.py`, `ea_all.py`/`.txt`, outputs `out.tar.gz`) |
| `searches/new10/` | the 4 new classes: LUTs, raw engine output, univariate forms (`univariate.py`), spectra, comparison with Helleseth-Ichanska-Kaleyski |
| `searches/fpgap/` | rerun of the search-only types with APN functions, printing every leaf, and EA check of every leaf: `search_job.sh`, `ea_job.sh`, `ea_leaves.py`, `run_all.sh`, `check_fp.py`, `status.txt`, `status.tar.xz` |
| `searches/t64/` | type 64 (`README.md`: exploratory, what it establishes): Hermitian criterion (`herm.py`), rank distributions (`rankdist64.py`), membership of the known classes (`iso64.py`), point-by-point solver (`t64.c`), switching (`switch64.py`, `bfs64.py`, `sw2.c`), sparse polynomials (`poly64.c`), structured search `tri10/` |
| `searches/t41/` | type 41: the split partial run (not needed: BBF, C and D prove it empty) |
| `lib/apnsearch/` | Python library: EA-equivalence with certificates (`equivalence.py`, options `kernels`, `transitive`), ortho-derivative invariants (`invariants.py`), C core |
| `data/quadratic_classes.tsv` | the 29 known quadratic APN classes of F_2^10 (25 published + 4 new): name, invariant, M, LUT |
| `DATA.md` | raw data not included, with sizes |
| `tools/check_manifest.py`, `MANIFEST.sha256` | SHA-256 of every tracked file |

## Building and running the searches

```sh
cd engine
cc -O3 -march=native -o sym10 sym10.c
cc -O3 -march=native -o sym10f73 sym10f73.c
../types/make_w.sh                      # W/w_ID.txt for all 299 types
./sym10 -dfs -dedupe < ../W/w_231.txt   # type 231: 961 leaves in a few seconds
```

The production runs used builds of `sym10.c` at successive stages under the names `sym10b` (`run_type.sh`, no
symmetry options), `sym10s` (`knuth_sym.sh`), `sym10cp` (`run_fast.sh`, `run_cfg.sh`, `run_strong.sh`,
`run_both.sh`), `sym10cp2` (`searches/t41/run41t.sh`), `sym10f` (type-47 Knuth estimates). The current
`sym10.c` accepts all their options; it reproduces the server node counts (type 231, the forced type-73 runs).
`sym10f73` was used for the type-73 rerun (`searches/t73/rerun_all/run73a.sh`) and the fingerprint-gap rerun.
The type-73 production run (192 pieces) used the options of `run_strong.sh` with `-shard r 192 -sdepth 10`
(see the header of `run73a.sh`).

Dedicated solvers (build lines are in the source headers):

```sh
cc -O3 -DM=5 -o full24g searches/t24/c/full24g.c          # type 24 (env GE, PART, MINKEY, OUT; see run24.sh)
cc -O3 -DKU=3 -DMZ=4 -o stage2i searches/t47/stage2i.c     # type 47 (env XONLY or PREFIX; see run47i.sh)
cc -O3 -DKU=3 -DMZ=4 -o stage2g searches/t47/stage2g.c
```

The run scripts were written for Linux servers. They wait while the load is high (`polite.sh` reads
`/proc/loadavg`) and run at nice 19. They expect the binaries and the input files in the working directory.
`searches/t24/run24.sh` was reconstructed from the job outputs; the other scripts are the originals.

## Fingerprint-gap rerun

The fingerprint-gap rerun (`searches/fpgap`) finished on 2026-11-08. `status.txt` and `status.tar.xz` hold its
final state: all 67 search jobs agree with the original runs, all 641 EA chunks ended, and all 4,851,099 leaves
are EA-equivalent to an apndb representative (0 failures). The leaf counts per type and class are in `status.txt`.

## Conventions

- Vectors of F_2^10 are integers; a function is a lookup table of 1,024 integers.
- Types are numbered as in `types/types10.tsv` (columns: id, p, kind, factors of A, factors of B, columns of A,
  columns of B).
- Classes are named by the ortho-derivative invariant (`types/classify.py`: SHA-1 of the differential and
  extended Walsh spectra of the ortho-derivative, first 16 hex digits). The invariant separates the 29 known
  classes. Equal invariants are turned into proved EA-equivalences by `lib/apnsearch/equivalence.py`, which
  returns (L1, L2); `check_certificate` then verifies that G + L1 o F o L2 is affine. In the stored runs this was done
  for every leaf (see mode (c)); `verify.py` recomputes it for a sample only.

## License

MIT (see `LICENSE`).
