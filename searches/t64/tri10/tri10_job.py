#!/usr/bin/env python3
"""One shard of the structured n = 10 search: run tri10, classify the printed functions (one per new cheap spectrum)
by the full ortho-derivative invariant, and test every class not seen before for CCZ-equivalence to a permutation.
Appends to results/classes.tsv (one line per class found in this shard: inv, known name or NEW, M, status,
config, lut) and to results/jobs.tsv (config, candidates, APN hits, spectra). Stop everything: touch STOP_TRI10.
usage: tri10_job.py BASE TERMS FIELD SHARD NSHARDS"""
import fcntl, hashlib, os, subprocess, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "lib"))
from apnsearch import ccz  # noqa: E402
from apnsearch.invariants import ortho_invariant  # noqa: E402

base, terms, field, shard, nsh = sys.argv[1:6]
root = here
if os.path.exists(os.path.join(root, "STOP_TRI10")):
    sys.exit(0)
out_dir = os.path.join(root, "results")
os.makedirs(out_dir, exist_ok=True)
cfg = f"base={base},terms={terms},field={field},shard={shard}/{nsh}"
r = subprocess.run([os.path.join(here, "tri10"), "-base", base, "-terms", terms, "-field", field, "-shard", shard,
                    "-nshards", nsh], capture_output=True, text=True)
stats = [l for l in r.stderr.splitlines() if l.startswith("done")]
known = {}
for line in open(os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")):
    t = line.split("\t")
    if t[0] != "name":
        known[t[1]] = t[0]
lines = [l for l in r.stdout.splitlines() if "\t" in l]
rows = []
for l in lines:
    F = np.array(list(map(int, l.split("\t")[1].split(","))))
    inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
    name = known.get(inv, "NEW")
    M = st = "-"
    if name == "NEW":
        st, info = ccz.ccz_to_permutation_dual(F, 10, quadratic=True)
        M = info.get("max_dim_NB")
    rows.append(f"{inv}\t{name}\t{M}\t{st}\t{cfg}\t{l.split(chr(9))[1] if name == 'NEW' else ''}\n")
    if st == "permutation":
        open(os.path.join(out_dir, "PERMUTATION"), "a").write(l + "\n")


def append(path, text):
    with open(path, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(text)
        fcntl.flock(f, fcntl.LOCK_UN)


append(os.path.join(out_dir, "classes.tsv"), "".join(rows))
append(os.path.join(out_dir, "jobs.tsv"), f"{cfg}\t{stats[0] if stats else 'no stats'}\tspectra={len(lines)}\n")
