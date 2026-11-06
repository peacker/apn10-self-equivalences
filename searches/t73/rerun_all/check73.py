#!/usr/bin/env python3
"""Check the type-73 rerun (shard J of 3072, J = r + 192 i) against the 2026-11-02 run (piece r of 192, t73_final/out5):
for each piece r with all 16 shards done: nodes = sum(nodes_J) - 15 * TOP (TOP = nodes above depth 10, counted in
every shard), leaves = sum(leaves_J), and the fingerprints printed by the old run are among the new ones."""
import glob, os, re, collections
here = os.path.dirname(os.path.abspath(__file__))
TOP = 7758789
new = {}
for f in glob.glob(os.path.join(here, "e?", "out", "a*.done")):
    j = int(os.path.basename(f)[1:-5]); m = re.search(r"nodes (\d+), leaves (\d+)", open(f).read())
    if m: new[j] = (int(m.group(1)), int(m.group(2)), os.path.splitext(f)[0] + ".out")
old = os.path.join(here, "..", "final", "out5")
ok = bad = pend = 0; nleaves = 0
for r in range(192):
    js = [r + 192 * i for i in range(16)]
    if not all(j in new for j in js): pend += 1; continue
    o = re.search(r"nodes (\d+), leaves (\d+), distinct fingerprints (\d+)", open(os.path.join(old, f"t_73.{r}.done")).read())
    on, ol, od = map(int, o.groups())
    nn = sum(new[j][0] for j in js) - 15 * TOP; nl = sum(new[j][1] for j in js)
    nk = collections.Counter(); 
    for j in js:
        for line in open(new[j][2]):
            if line.startswith("K"): nk[line.split()[0]] += 1
    ok_k = set(l.split()[0] for l in open(os.path.join(old, f"t_73.{r}.out")) if l.startswith("K")) == set(nk)
    good = nn == on and nl == ol and sum(nk.values()) == nl and ok_k and len(nk) == od
    ok += good; bad += not good; nleaves += nl
    if not good: print(f"piece {r}: old nodes {on} leaves {ol} fp {od}; new nodes {nn} leaves {nl} printed {sum(nk.values())} fp {len(nk)} same-fp-set {ok_k}")
print(f"pieces checked {ok + bad}: {ok} agree, {bad} differ; {pend} not complete; leaves in complete pieces {nleaves}")
