#!/usr/bin/env python3
"""Type 24: check the stored outputs of the full24g run (results/t24_out.tar.gz, unpacked to results/out/).
Every one of the 4,096 jobs (Gold exponent GE = 3, 5 and part R = 0..2047 of 2048) must have a .done file with
"rc 0" and a summary line in its .out; the leaves are the lines of the .lut files ("tag LUT" or a bare LUT).
Prints the number of jobs, the summed run time, the leaves, and for each leaf the APN test and the
ortho-derivative invariant with its apndb name. Expected: 4,096 jobs, 3 leaves (g3.254 Gold x^3, g3.558 F4 a=1,
g5.1838 Gold x^9). The EA certificates are in results/ea_check24.py / ea_check24.txt."""
import glob, hashlib, os, re, sys, tarfile
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
res = os.path.join(here, "results")
out = os.path.join(res, "out")
if not os.path.isdir(out):
    with tarfile.open(os.path.join(res, "t24_out.tar.gz")) as t:
        try:
            t.extractall(res, filter="data")
        except TypeError:
            t.extractall(res)
known = {t[1]: t[0] for t in (l.rstrip("\n").split("\t") for l in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")))}


def is_apn(F):
    for a in range(1, 1024):
        d = F ^ F[np.arange(1024) ^ a]
        if np.bincount(d, minlength=1024).max() > 2:
            return False
    return True


want = {f"g{ge}.{r}" for ge in (3, 5) for r in range(2048)}
ok_jobs, secs, bad = set(), 0, []
for d in glob.glob(os.path.join(out, "*.done")):
    tag = os.path.basename(d)[:-5]; s = open(d).read()
    o = open(d[:-5] + ".out").read()
    if s.startswith("rc 0") and re.search(r"leaves \d+ APN \d+", o):
        ok_jobs.add(tag); secs += int(re.search(r"secs (\d+)", s).group(1))
    else:
        bad.append(tag)
leaves = []
for f in sorted(glob.glob(os.path.join(out, "*.lut"))):
    for line in open(f):
        w = line.split()
        if not w: continue
        lut = w[1:] if len(w) > 2 else (w[1] if len(w) == 2 else w[0]).split(",")
        leaves.append((os.path.basename(f)[:-4], np.array(list(map(int, lut)))))
print(f"jobs finished with rc 0 and a summary: {len(ok_jobs & want)} of {len(want)} (bad {len(bad)}); "
      f"summed run time {secs / 3600:.0f} core-hours; leaves {len(leaves)}")
names = []
for tag, F in leaves:
    inv = hashlib.sha1(repr(ortho_invariant(F, 10)).encode()).hexdigest()[:16]
    apn = is_apn(F)
    names.append((tag, known.get(inv, "NOT IN APNDB"), apn))
    print(f"  {tag}: APN {apn}, invariant {inv} = {known.get(inv, 'NOT IN APNDB')}")
ok = ok_jobs == want and sorted(names) == [("g3.254", "Gold x^3", True), ("g3.558", "F4 a=1", True), ("g5.1838", "Gold x^9", True)]
print("TYPE 24 RUN OK" if ok else "TYPE 24 RUN FAILED")
sys.exit(0 if ok else 1)
