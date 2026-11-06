#!/usr/bin/env python3
"""Class-level BFS over codimension-1 Hermitian switching (switch64.neighbours) from the 16 known classes of
type 64 (seeds64.txt), in parallel. Expands up to REPS representatives per class. New classes go to OUT.new_luts.
usage: bfs64.py [--reps K] [--procs P] [--rounds R] OUT_PREFIX"""
import os, sys
from multiprocessing import Pool
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
from switch64 import neighbours, lut_of, inv_hash

def expand(V):
    res = {}
    for W in neighbours(V, max_new=8):
        F = lut_of(W); h = inv_hash(F)
        if h not in res: res[h] = (W, F.tolist())
    return res

if __name__ == "__main__":
    a = sys.argv[1:]; reps, procs, rounds = 2, 8, 10
    for opt in ("--reps", "--procs", "--rounds"):
        if opt in a:
            i = a.index(opt); v = int(a[i + 1]); del a[i:i + 2]
            if opt == "--reps": reps = v
            elif opt == "--procs": procs = v
            else: rounds = v
    out = a[0]
    known = {}
    for line in open(os.path.join(here, "..", "..", "data", "quadratic_classes.tsv")):
        t = line.rstrip("\n").split("\t")
        if t[0] != "name": known[t[1]] = t[0]
    classes = {}
    for line in open(os.path.join(here, "seeds64.txt")):
        name, codes = line.rstrip("\n").split("\t"); V = list(map(int, codes.split()))
        classes.setdefault(inv_hash(lut_of(V)), []).append(V)
    expanded = set()
    log = open(out + ".log", "a")
    with Pool(procs) as pool:
        for rnd in range(rounds):
            todo = [(h, V) for h, vs in classes.items() for V in vs[:reps] if (h, tuple(V)) not in expanded]
            if not todo: break
            for h, V in todo: expanded.add((h, tuple(V)))
            nnew = 0
            for res in pool.imap_unordered(expand, [V for _, V in todo]):
                for h, (W, F) in res.items():
                    if h not in classes:
                        classes[h] = [W]; nnew += 1
                        name = known.get(h, "NEW")
                        msg = "round %d class %s %s V=%s" % (rnd, h, name, " ".join(map(str, W)))
                        print(msg, flush=True); log.write(msg + "\n"); log.flush()
                        if name == "NEW":
                            with open(out + ".new_luts", "a") as fo: fo.write(",".join(map(str, F)) + "\n")
                    elif len(classes[h]) < reps:
                        classes[h].append(W)
            print("round %d: %d classes, %d new, %d known-named" % (rnd, len(classes), nnew,
                  sum(1 for h in classes if h in known)), flush=True)
