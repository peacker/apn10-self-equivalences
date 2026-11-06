"""Compare the new EA options (kernels=True; transitive=True for Gold classes) with ea_all.txt on type-73 leaves of
every class: the result must again be 'equivalent' with a checked certificate. Prints time per class."""
import collections, glob, hashlib, os, sys, time
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "lib"))
from apnsearch.equivalence import ea_equivalence, check_certificate  # noqa: E402
tsv = os.path.join(here, "..", "..", "..", "data", "quadratic_classes.tsv")
known = {t[1]: t for t in (l.rstrip("\n").split("\t") for l in open(tsv)) if t[0] != "name"}
done = {l.split()[0]: l.split()[2] for l in open(os.path.join(here, "ea_all.txt")) if ": equivalent" in l}
per = int(sys.argv[1]) if len(sys.argv) > 1 else 2
cnt = collections.Counter(); tm = collections.defaultdict(float); bad = 0
for f in sorted(glob.glob(os.path.join(here, "e?", "out", "a*.out"))):
    for line in open(f):
        if not line.startswith("K"): continue
        s = line.split()[1]; h = hashlib.sha1(s.encode()).hexdigest()[:16]
        inv = done.get(h)
        if inv is None or cnt[inv] >= per: continue
        F = np.array(list(map(int, s.split(",")))); D = np.array(list(map(int, known[inv][-1].split(","))))
        g = known[inv][0].startswith("Gold")
        t0 = time.time(); c = ea_equivalence(D, F, 10, transitive=g, kernels=True)
        ok = c is not None and check_certificate(D, F, 10, *c); tm[inv] += time.time() - t0; cnt[inv] += 1
        bad += not ok
for inv in cnt: print(f"{known[inv][0]:45s} {cnt[inv]} leaves, {tm[inv] / cnt[inv]:.1f}s each")
print("all equivalent with checked certificate" if not bad else f"{bad} FAILED")
