"""Time the EA test with and without transitive=True on Gold leaves of the type-73 rerun (already checked by ea_all.py)."""
import glob, hashlib, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "lib"))
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402
tsv = os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")
known = {t[1]: t for t in (l.rstrip("\n").split("\t") for l in open(tsv)) if t[0] != "name"}
gold = {l.split()[0]: l.split()[2] for l in open(os.path.join(here, "ea_all.txt")) if "Gold" in l}
want = int(sys.argv[1]) if len(sys.argv) > 1 else 6
old = "--old" in sys.argv
got = 0
for f in sorted(glob.glob(os.path.join(here, "e?", "out", "a*.out"))):
    for line in open(f):
        if not line.startswith("K"): continue
        s = line.split()[1]; h = hashlib.sha1(s.encode()).hexdigest()[:16]
        if h not in gold: continue
        F = np.array(list(map(int, s.split(",")))); D = np.array(list(map(int, known[gold[h]][-1].split(","))))
        t0 = time.time(); c = ea_equivalence(D, F, 10, transitive=True, kernels=True); t1 = time.time()
        ok = c is not None and check_certificate(D, F, 10, *c)
        msg = f"{h} {known[gold[h]][0]} transitive: {ok} {t1 - t0:.1f}s"
        if old:
            t0 = time.time(); c2 = ea_equivalence(D, F, 10); msg += f"; old: {c2 is not None} {time.time() - t0:.0f}s"
        print(msg, flush=True); got += 1
        if got >= want: sys.exit()
