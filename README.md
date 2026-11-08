# Quadratic APN functions of F_2^10 with a prime-order self-equivalence: code and data

Code, generated data and stored results of the paper

> E. Bellini, *Quadratic APN functions of F_2^10 with a prime-order self-equivalence* (title provisional),
> in preparation for IEEE Transactions on Information Theory.

The paper classifies the quadratic APN functions F of F_2^10 that have a linear self-equivalence (A, B) of prime
order, i.e. beta(Ax, Ay) = B beta(x, y) for the polar form beta of F. Up to conjugation and powers there are 299
such types (A, B). Each type is settled by a proof criterion, by an exhaustive search, or by both. Type 64 stays
open. The repository contains everything needed to re-check each claim. `verify.py` re-runs the cheap parts and
checks the stored outputs of the expensive searches in about 5 minutes. The searches themselves took about
40,000 core-hours and are not re-run by `verify.py`.

## Requirements

- Python 3.10 or later with NumPy (`pip install -r requirements.txt`)
- a C compiler (`cc`; the C core of `lib/apnsearch` is compiled on first import into `~/.cache/apnsearch`,
  override with `APNSEARCH_BUILD=<dir>`)

Tested on macOS 26 (Apple clang 16, Python 3.13, NumPy 2.2). The searches were run on Linux servers
(gcc, `nice 19`, see `engine/polite.sh`).

## Quick start

```sh
python3 verify.py              # all quick checks, about 5 minutes on one core
python3 verify.py --skip-ea    # without recomputing EA certificates
```

`verify.py` prints one line per check and `ALL CHECKS PASSED` at the end (exit status 1 otherwise). It unpacks
the `.tar.gz` archives of stored outputs next to themselves on first use (git ignores the unpacked copies).
Last run: all 43 checks passed, 270 s from a fresh clone (macOS, Apple M5 Pro, one core).

## Claims and how to re-check them

| Claim (paper) | Files | How to re-check | Cost |
|---|---|---|---|
| There are 299 types (A, B) of prime order at n = 10 (60 at n = 6, 107 at n = 8) | `types/enumerate.py`, `types/types{6,8,10}.tsv` | `python3 types/enumerate.py 10` (in `verify.py`) | seconds |
| The space W of polar forms of each type | `types/sym_setup.py`, `types/make_w.sh`; `engine/W3/` (the 35 spaces used by the searches) | `types/make_w.sh` writes `W/w_ID.txt` for all 299 types; deterministic, equal to `engine/W3` (in `verify.py`) | 1 min |
| Proof criteria C (compression), R (restriction + Chevalley-Warning), Q (order 3), W=0 prove 276 types | `proofs/proofs10.py` -> `proofs/proofs10.txt` | `python3 proofs/proofs10.py W` (in `verify.py`) | 10 s |
| Wong criterion proves types 130, 169 | `proofs/wong10.py` -> `proofs/wong10.txt` | `python3 proofs/wong10.py W 130 169` (in `verify.py`) | 1 s |
| Odd-dimension lemma ("lemma56") excludes types 55-59, in particular 55, 56 | `proofs/lemma56/check56.py` (numerical check of the lemma) | `python3 proofs/lemma56/check56.py 56` (in `verify.py`) | 1 s |
| Dual intertwiner lemma D excludes 95 types | `proofs/lemma_dual.py` -> `proofs/lemma_dual.txt` | in `verify.py` | 1 s |
| The bound of Baudrin et al. (DCC 2025, Cor. 3) excludes the 65 types with dim Fix(A) + dim Fix(B) >= 12 (among them type 41) | `proofs/bbf.py` -> `proofs/bbf_excluded.txt` | in `verify.py` | 1 s |
| Master table: 242 types by proof and search, 38 by proof only, 18 by search only, type 64 open | `proofs/master_table.py` (inputs: `types/types10.tsv`, `proofs10.txt`, `lemma_dual.txt`, `bbf_excluded.txt`, `search_status.txt`, `server_done.txt`) -> `proofs/master_table.tsv` | in `verify.py` | 1 s |
| The generic search engine is exhaustive and reproducible | `engine/sym10.c`, `engine/sym10f73.c`, `engine/validation/` | type 231 re-run (91,636 nodes, 961 leaves, as on the server); 12 forced type-73 runs (in `verify.py`) | 20 s |
| Search-only types with APN functions (129, 170, 171, 175, 176, 180, 231): every leaf printed again and EA-tested | `searches/fpgap/` | `python3 searches/fpgap/check_fp.py` (in `verify.py`); full re-run: `run_all.sh` | rerun complete (2026-11-08): 67 of 67 jobs agree, all 4,851,099 leaves EA-equivalent to a known class; search 9 core-hours, EA checks 1,332 core-hours; check 2 s |
| **Type 24** (coefficients in F_32): exactly 3 classes, Gold x^3, x^9 and x^3 + Tr(x^9) | `searches/t24/` (normal form `form24.py`, plane lemma, solver `c/full24g.c`, `run24.sh`, outputs `results/t24_out.tar.gz`, EA check `results/ea_check24.txt`) | `python3 searches/t24/check_run.py`; EA certificates recomputed in `verify.py`; full re-run: 4,096 jobs of `run24.sh` | 7,500 core-hours; check 1 min |
| **Type 47**: empty | `searches/t47/` (stage 1 `stage1.py`, `stage1b.py` -> 806 classes `cls/`; rank-h lemma; solvers `stage2g.c`, `stage2i.c`; outputs `results/t47_out.tar.gz`) | `python3 searches/t47/check_coverage.py` (in `verify.py`): 703 classes dead, 8,056 slices of the other 103 all searched, no leaf | 13,261 core-hours; check 1 min |
| **Type 73** (coefficients in F_4): exactly 24 classes, 20 published and 4 new | `searches/t73/final/` (first run, 192 pieces), `searches/t73/rerun_all/` (rerun printing every leaf, 3,072 shards; `ea_all.txt`: every one of the 4,536 leaves EA-tested) | `check73.py` (rerun = first run, piece by piece) and the counts of `ea_all.txt` in `verify.py`; one EA certificate per class recomputed; full EA re-check: `python3 searches/t73/rerun_all/ea_all.py` (about 30 min) | 12,600 + 5,800 (rerun) core-hours |
| The 4 new classes are APN, quadratic, and differ from all 25 published classes | `data/quadratic_classes.tsv` (copy of apndb `n10/quadratic_classes.tsv`), `searches/new10/` (univariate forms over F_4; shortest forms in `short_forms.txt`, checked by `check_short_forms.py`) | in `verify.py` (APN test, degree, ortho-derivative spectra) | seconds |
| Type 41 is empty | BBF (`proofs/bbf.py`), criteria C and D; the partial search is in `searches/t41/` (not needed) | in `verify.py` | - |
| **Type 64** (canonical triplicates) is open; 16 of the 29 known classes lie in it; four searches find nothing else | `searches/t64/` (`iso64.py`, `herm.py`, `t64.c`, `switch64.py`, `bfs64.py`, `sw2.c`, `poly64.c`, `tri10/`) | `python3 searches/t64/iso64.py` (in `verify.py`); `python3 searches/t64/herm.py` (Hermitian criterion on x^3, x^9, x^3 + Tr(x^9)) | 1 min |

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
| `types/` | `enumerate.py` (types of prime-order (A, B), Beierle-Leander style), `types{6,8,10}.tsv`, `sym_setup.py` (the space W and the centraliser bases, input of `sym10`), `make_w.sh`, `reorder.py`, `classify.py` (name LUTs by the ortho-derivative invariant) |
| `proofs/` | `proofs10.py`/`.txt`, `wong10.py`/`.txt`, `lemma_dual.py`/`.txt`, `lemma56/check56.py`, `bbf.py`/`bbf_excluded.txt`, `master_table.py`/`.tsv`, `search_status.txt` (per type: which server runs are complete), `server_done.txt` (first line of every server `.done` file; machine names replaced by `srv1`, `srv2`) |
| `engine/` | `sym10.c` (generic symmetric DFS over W; usage in the header), `sym10f73.c` (same search, faster checkpoints, `-all` prints every leaf), run scripts (`run_type.sh`, `run_fast.sh`, `run_cfg.sh`, `run_strong.sh`, `run_both.sh`, `polite.sh`), Knuth estimates (`costs10.tsv`, `knuth_both.tsv`, `knuth_sym.tsv`, ...), `W3/` (W spaces of the hard types), batch files, `intertwine.py`, `gen_force.py`, `force_*.txt`, `validation/` (runs on known types, n = 6, 8 and forced type-73 prefixes) |
| `searches/t24/` | type 24: `form24.py`, `topbent.py`, `sparse_sub.py`, `c/full24.c` (n = 6 validation), `c/full24g.c` (production), `run24.sh`, `check_run.py`, `results/` |
| `searches/t47/` | type 47: `stage1.py`, `stage1b.py`, `cls/` (806 classes), `s2/` (rank-h filter), `stage2.c`, `stage2g.c`, `stage2i.c`, `gold6.py` (n = 6 validation), job lists, run scripts, `check_coverage.py`, `results/` |
| `searches/t73/` | type 73: `final/` (first complete run: `classes.txt`, `ea_check.py`/`.txt`, LUTs, outputs `out5.tar.gz`), `rerun_all/` (rerun with every leaf: `run73a.sh`, `check73.py`, `ea_all.py`/`.txt`, outputs `out.tar.gz`) |
| `searches/new10/` | the 4 new classes: LUTs, raw engine output, univariate forms (`univariate.py`), spectra, comparison with Helleseth-Ichanska-Kaleyski |
| `searches/fpgap/` | rerun of the search-only types with APN functions, printing every leaf, and EA check of every leaf: `search_job.sh`, `ea_job.sh`, `ea_leaves.py`, `run_all.sh`, `check_fp.py`, `status.txt`, `status.tar.xz` |
| `searches/t64/` | type 64: Hermitian criterion (`herm.py`), membership of the known classes (`iso64.py`), point-by-point solver (`t64.c`), switching (`switch64.py`, `bfs64.py`, `sw2.c`), sparse polynomials (`poly64.c`), structured search `tri10/` |
| `searches/t41/` | type 41: the split partial run (not needed: BBF, C and D prove it empty) |
| `lib/apnsearch/` | Python library: EA-equivalence with certificates (`equivalence.py`, options `kernels`, `transitive`), ortho-derivative invariants (`invariants.py`), C core |
| `data/quadratic_classes.tsv` | the 29 known quadratic APN classes of F_2^10 (25 published + 4 new): name, invariant, M, LUT |
| `DATA.md` | raw data not included, with sizes |

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
  returns (L1, L2) and checks the certificate independently.

## License

MIT (see `LICENSE`).
