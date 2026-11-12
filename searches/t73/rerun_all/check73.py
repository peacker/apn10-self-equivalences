#!/usr/bin/env python3
"""Check the type-73 rerun (shard J of 3072, J = r + 192 i) against the 2026-11-02 run (piece r of 192, final/out5).

Exact identity checks (duplicates never compensate for missing records):
 - the first run has exactly the pieces r = 0..191 (final/out5/t_73.r.done), each with a "done:" line, no node-limit
   marker, and the progress file final/out5/t_73.r.prog ending with DONE;
 - jobs_all.txt lists exactly the shards 0..3071; the rerun (e1/out, e4/out) has exactly one record aJ.done for each
   of them (no shard missing, none duplicated), each "rc 0" with a "done:" line and no node-limit marker;
 - for each piece r: nodes = sum(nodes_J) - 15 * TOP (TOP = nodes above depth 10, counted in every shard),
   leaves = sum(leaves_J), and the fingerprints printed by the first run are those of the rerun;
 - the leaves printed by the rerun ("K<fp> <lut>" lines of aJ.out) are pairwise distinct, and the set of their
   hashes sha1(lut)[:16] is exactly the set of leaf hashes in ea_all.txt (one line per leaf, no duplicates), all
   "equivalent".
Prints a summary; exit status 1 on any failure."""
import collections, glob, hashlib, os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
TOP = 7758789
NPIECE, NSHARD = 192, 3072
errors = []
LIMIT_MARK = ("node limit", "INCOMPLETE")

# first run
old = os.path.join(here, "..", "final", "out5")
old_ids = sorted(int(re.match(r"t_73\.(\d+)\.done$", os.path.basename(f)).group(1))
                 for f in glob.glob(os.path.join(old, "t_73.*.done")))
if old_ids != list(range(NPIECE)):
    errors.append(f"first run: piece ids are not exactly 0..{NPIECE - 1} ({len(old_ids)} records)")
for r in range(NPIECE):
    p = os.path.join(old, f"t_73.{r}.done")
    s = open(p).read() if os.path.exists(p) else ""
    pr = os.path.join(old, f"t_73.{r}.prog")
    prog = open(pr).read() if os.path.exists(pr) else ""
    if not re.search(r"done: nodes \d+, leaves \d+, distinct fingerprints \d+", s) or any(m in s for m in LIMIT_MARK) \
            or not prog.rstrip().endswith("DONE"):
        errors.append(f"first run piece {r}: incomplete record")

# rerun
jobs = sorted(int(l) for l in open(os.path.join(here, "jobs_all.txt")) if l.strip())
if jobs != list(range(NSHARD)):
    errors.append("jobs_all.txt is not exactly 0..3071")
recs = collections.defaultdict(list)
for f in glob.glob(os.path.join(here, "e?", "out", "a*.done")):
    recs[int(os.path.basename(f)[1:-5])].append(f)
if sorted(recs) != list(range(NSHARD)):
    errors.append(f"rerun: shard ids are not exactly 0..{NSHARD - 1} ({len(recs)} distinct ids)")
dups = [j for j, v in recs.items() if len(v) > 1]
if dups:
    errors.append(f"rerun: {len(dups)} shards have more than one record, e.g. {dups[:5]}")
new = {}
for j, fs in recs.items():
    s = open(fs[0]).read()
    m = re.match(r"rc (\d+) done: nodes (\d+), leaves (\d+)", s)
    if not m or m.group(1) != "0" or any(x in s for x in LIMIT_MARK):
        errors.append(f"rerun shard {j}: not complete ({s.strip()[:80]})"); continue
    new[j] = (int(m.group(2)), int(m.group(3)), fs[0][:-5] + ".out")

ok = bad = pend = 0; nleaves = 0
emitted = collections.Counter()
for r in range(NPIECE):
    js = [r + NPIECE * i for i in range(16)]
    if not all(j in new for j in js):
        pend += 1; continue
    o = re.search(r"nodes (\d+), leaves (\d+), distinct fingerprints (\d+)", open(os.path.join(old, f"t_73.{r}.done")).read())
    on, ol, od = map(int, o.groups())
    nn = sum(new[j][0] for j in js) - 15 * TOP; nl = sum(new[j][1] for j in js)
    nk = collections.Counter()
    for j in js:
        for line in open(new[j][2]):
            if line.startswith("K"):
                fp, lut = line.split()[:2]
                nk[fp] += 1
                emitted[hashlib.sha1(lut.encode()).hexdigest()[:16]] += 1
    ok_k = set(l.split()[0] for l in open(os.path.join(old, f"t_73.{r}.out")) if l.startswith("K")) == set(nk)
    good = nn == on and nl == ol and sum(nk.values()) == nl and ok_k and len(nk) == od
    ok += good; bad += not good; nleaves += nl
    if not good:
        print(f"piece {r}: old nodes {on} leaves {ol} fp {od}; new nodes {nn} leaves {nl} printed {sum(nk.values())} "
              f"fp {len(nk)} same-fp-set {ok_k}")
print(f"pieces checked {ok + bad}: {ok} agree, {bad} differ; {pend} not complete; leaves in complete pieces {nleaves}")
if bad or pend:
    errors.append(f"{bad} pieces differ, {pend} not complete")

# EA records <-> emitted leaves
ea = [l.split() for l in open(os.path.join(here, "ea_all.txt")) if l.strip()]
ea_h = collections.Counter(x[0] for x in ea)
rep = [h for h, c in emitted.items() if c > 1]
if rep:
    errors.append(f"{len(rep)} leaf LUTs emitted more than once")
if [h for h, c in ea_h.items() if c > 1]:
    errors.append("ea_all.txt has duplicate leaf records")
if set(ea_h) != set(emitted):
    errors.append(f"ea_all.txt leaf set differs from the emitted leaves ({len(set(emitted) - set(ea_h))} missing, "
                  f"{len(set(ea_h) - set(emitted))} extra)")
if not all(" ".join(x).count(": equivalent") == 1 for x in ea):
    errors.append("ea_all.txt has a leaf that is not 'equivalent'")
print(f"identity: first-run pieces {len(old_ids)}, rerun shards {len(recs)} (duplicates {len(dups)}), emitted leaves "
      f"{sum(emitted.values())} ({len(emitted)} distinct), EA records {len(ea)}, EA set = emitted set: {set(ea_h) == set(emitted)}")
for e in errors:
    print("  !", e)
print("TYPE 73 CHECK OK" if not errors else "TYPE 73 CHECK FAILED")
sys.exit(1 if errors else 0)
