# Data not included

The repository keeps the outputs needed to check the claims. Large raw outputs are stored compressed, and
material that no claim of the paper depends on is left out. Paths below are relative to the author's working
tree `learning/sym10/` (and `learning/quad/`), unless a server is named. Server outputs are only on the servers
and in the author's archive. They are available on request.

## Included in compressed form

| Archive in this repository | Unpacked size | Content |
|---|---|---|
| `searches/fpgap/status.tar.xz` (45 MB) | 328 MB | `status/out` (67 search `.done` files) and `status/ea` (one line per printed leaf: hash, invariant, class, result) of the fingerprint-gap rerun, final state (complete on 2026-11-08) |
| `searches/t73/rerun_all/out.tar.gz` (3.4 MB) | 42 MB | the 3,072 shard outputs of the type-73 rerun (every leaf as a LUT) |
| `searches/t73/final/out5.tar.gz` (1 MB) | 4.3 MB | the 192 piece outputs of the first complete type-73 run |
| `searches/t47/results/t47_out.tar.gz` (1.3 MB) | 49 MB | outputs of all finished type-47 jobs (stage2g `out2`, stage2i `outi`, synced from two servers) |
| `searches/t24/results/t24_out.tar.gz` (1.2 MB) | 32 MB on disk (0.7 MB of text) | the 4,096 job outputs of the type-24 run (`.out`, `.done`, `.lut`, `.err`) |
| `searches/t41/out2_partial.tar.gz` | 0.8 MB | partial type-41 split run (not needed for any claim) |

## Not included

| Data | Size | Where it lives | Why not included |
|---|---|---|---|
| Printed leaves of the fingerprint-gap rerun (`out/t_*.out.gz`, one LUT per leaf, 4,851,099 leaves) | 7.8 GB (201 files, checked on the server on 2026-11-13) | one server (`fpgap/out/`); available on request | too large; every leaf is EA-checked in `status/ea` (included); `search_job.sh` regenerates them (about 9 core-hours), and re-checking all equivalences then takes about 1,300 core-hours |
| Raw `.out`/`.err` files of the generic `sym10` runs over the 299 types | not copied locally | servers (`sym10/out`, `sym10b/out2` ... `out5`) | summarized in `proofs/server_done.txt` (first line of every `.done` file: nodes, leaves, seconds) and `proofs/search_status.txt`; the leaves of the types with APN functions are re-printed by `searches/fpgap` |
| First type-73 run (32 pieces, `out4/old73`, 381 core-hours) | not copied locally | server | superseded by the 192-piece run and the rerun (included) |
| `t24/c/q75/` | 369 MB | `learning/sym10/t24/c/q75` | Kim-like permutation side study (75 affine classes of quadratic permutations of F_2^5); not used by the type-24 classification |
| `t24/c/L_*.txt`, `walk*.txt`, `kl/`, `vc/` and the other Kim-like programs (`kimlike.c`, `vcheck.c`, `samp24.c`, ...) | 46 MB | `learning/sym10/t24/c/` | same side study; not used by the classification |
| Earlier type-24 solvers `full24b.c` ... `full24f.c` | < 0.2 MB | `learning/sym10/t24/c/` | superseded by `full24g.c` (included with `full24.c`) |
| `t47/clsr/`, `t47/tim*/`, `t47/r3/` (except `stab.c`, `covered_u.txt`), `t47/out2_e1/`, profiling variants of the solver | 5 MB | `learning/sym10/t47/` | scratch: alternative class files, timing runs, profiling builds |
| `P73/` | 0.6 MB | `learning/sym10/P73/` | pencil decomposition of type 73, tried and not used |
| `W2/`, `W5/`, `W8/` | 1 MB | `learning/sym10/` | W spaces of earlier engine versions or of n = 8; regenerable with `types/sym_setup.py` (`W3/` is included) |
| `out5/` (top level) | 0.3 MB | `learning/sym10/out5/` | partial local copy of `t73_final/out5` |
| Compiled binaries (`sym10cp`, `sym10f`, `sym10k`, `sym10s`, solver builds) | < 1 MB | `learning/sym10/` | rebuilt from the sources |
| `sym10/fpgap/lib`, `fpgap/quadratic_classes.tsv` | 0.2 MB | `learning/sym10/fpgap/` | identical copies of `lib/apnsearch` and `data/quadratic_classes.tsv` |
| Type-64 structured search: the lookup tables printed per shard | - | server | only the classified results (`searches/t64/tri10/results/`) were kept |
