#!/usr/bin/env python3
"""Compare the rerun (status/out/*.done) with the original runs (../../proofs/server_done.txt: nodes, leaves per job) and
summarize the EA results (status/ea/*.txt). Prints one summary line, and every disagreement or failure."""
import glob, os, re, collections
here = os.path.dirname(os.path.abspath(__file__))
orig = {}
for l in open(os.path.join(here, "..", "..", "proofs", "server_done.txt")):
    m = re.match(r"(\S+) sym10/out/(t_\d+(?:\.\d+)?)\.done done: nodes (\d+), leaves (\d+)", l)
    if m and (m.group(1) == "srv1" or m.group(2) in ("t_170", "t_171", "t_231")):
        orig[m.group(2)] = (int(m.group(3)), int(m.group(4)))
agree = differ = 0; bad = []
for f in glob.glob(os.path.join(here, "status", "out", "*.done")):
    tag = os.path.basename(f)[:-5]; s = open(f).read()
    m = re.search(r"rc (\d+) done: nodes (\d+), leaves (\d+)", s)
    if not m or m.group(1) != "0": bad.append(f"{tag}: {s.strip()}"); continue
    if orig.get(tag) == (int(m.group(2)), int(m.group(3))): agree += 1
    else: differ += 1; bad.append(f"{tag}: new {m.group(2)} {m.group(3)} vs old {orig.get(tag)}")
res = collections.Counter(); classes = collections.defaultdict(collections.Counter); ends = 0
for f in glob.glob(os.path.join(here, "status", "ea", "*.txt")):
    t = os.path.basename(f).split(".")[0][2:]
    for l in open(f):
        if l.startswith("END"): ends += 1; continue
        if l.startswith("L"): continue
        r = "equivalent" if ": equivalent" in l else "FAIL"
        res[r] += 1
        if r == "FAIL": bad.append(l.strip())
        else: classes[t][l.split(": ")[0].split(" ", 2)[2]] += 1
nd = len(glob.glob(os.path.join(here, "status", "ea", "*.done")))
print(f"search jobs: {agree} agree, {differ} differ (of 67); EA chunks done {nd}, ended {ends}; leaves equivalent "
      f"{res['equivalent']}, failed {res['FAIL']}")
for t in sorted(classes): print(f"  type {t}: " + "; ".join(f"{c} {n}" for c, n in sorted(classes[t].items())))
for b in bad[:20]: print("  !", b)
