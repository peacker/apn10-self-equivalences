#!/usr/bin/env python3
"""Check the fingerprint-gap rerun (status/, unpacked from status.tar.xz) with exact identity sets.

 - search jobs: exactly the 67 tags of run_all.sh (t_T.S for T in 129, 175, 176, 180 and S = 0..15; t_170, t_171,
   t_231), one record each, "rc 0" with a "done:" line and no node-limit marker; nodes and leaves equal to the
   original run (../../proofs/server_done.txt, srv1 runs; also without a node-limit marker);
 - EA chunks: for each tag exactly the chunk ids 0..K-1 of run_all.sh (K = 12 for 175, 176, 180; 24 for 170, 171;
   1 otherwise), 641 in total, each "rc 0"; chunk c of K takes the leaves k = c, c + K, ... of the tag
   (ea_leaves.py), so it must end with "END L" (L = the leaves of the search job) and hold exactly
   len(range(c, L, K)) leaf records, all "equivalent", with pairwise distinct leaf hashes within the tag.
Prints one summary line, the classes per type, and every failure; exit status 1 on any failure."""
import collections, glob, os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
st = os.path.join(here, "status")
LIMIT_MARK = ("node limit", "INCOMPLETE")
TAGS = [f"t_{t}.{s}" for t in (129, 175, 176, 180) for s in range(16)] + ["t_170", "t_171", "t_231"]
K = lambda tag: 12 if tag.startswith(("t_175.", "t_176.", "t_180.")) else 24 if tag in ("t_170", "t_171") else 1
bad = []

orig = {}
for l in open(os.path.join(here, "..", "..", "proofs", "server_done.txt")):
    m = re.match(r"(\S+) sym10/out/(t_\d+(?:\.\d+)?)\.done done: nodes (\d+), leaves (\d+)", l)
    if m and (m.group(1) == "srv1" or m.group(2) in ("t_170", "t_171", "t_231")) and m.group(2) in TAGS:
        if any(x in l for x in LIMIT_MARK):
            bad.append(f"original {m.group(2)}: node limit marker")
        orig[m.group(2)] = (int(m.group(3)), int(m.group(4)))

recs = collections.Counter(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(st, "out", "*.done")))
if set(recs) != set(TAGS):
    bad.append(f"search records are not exactly the 67 tags (missing {sorted(set(TAGS) - set(recs))[:5]}, "
               f"extra {sorted(set(recs) - set(TAGS))[:5]})")
agree = differ = 0; leaves = {}
for tag in TAGS:
    p = os.path.join(st, "out", tag + ".done")
    if not os.path.exists(p):
        differ += 1; continue
    s = open(p).read()
    m = re.search(r"rc (\d+) done: nodes (\d+), leaves (\d+)", s)
    if not m or m.group(1) != "0" or any(x in s for x in LIMIT_MARK):
        bad.append(f"{tag}: {s.strip()}"); differ += 1; continue
    leaves[tag] = int(m.group(3))
    if orig.get(tag) == (int(m.group(2)), int(m.group(3))):
        agree += 1
    else:
        differ += 1; bad.append(f"{tag}: new {m.group(2)} {m.group(3)} vs old {orig.get(tag)}")

want_chunks = {f"{t}.{c}" for t in TAGS for c in range(K(t))}
have_done = {os.path.basename(f)[:-5] for f in glob.glob(os.path.join(st, "ea", "*.done"))}
have_txt = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(st, "ea", "*.txt"))}
if have_done != want_chunks or have_txt != want_chunks:
    bad.append(f"EA chunk ids are not exactly the {len(want_chunks)} expected ones "
               f"(missing {sorted(want_chunks - (have_done & have_txt))[:5]}, extra {sorted((have_done | have_txt) - want_chunks)[:5]})")
res = collections.Counter(); classes = collections.defaultdict(collections.Counter); ends = 0
hashes = collections.defaultdict(collections.Counter)
for ch in sorted(want_chunks & have_done & have_txt):
    tag, c = ch.rsplit(".", 1); c = int(c); k = K(tag)
    d = open(os.path.join(st, "ea", ch + ".done")).read()
    if not d.startswith("rc 0"):
        bad.append(f"EA chunk {ch}: {d.strip()}")
    n = 0; end = None
    for l in open(os.path.join(st, "ea", ch + ".txt")):
        if l.startswith("END"):
            end = int(l.split()[1]); ends += 1; continue
        if l.startswith("L"):
            continue
        n += 1
        r = "equivalent" if ": equivalent" in l else "FAIL"
        res[r] += 1
        hashes[tag][l.split()[0]] += 1
        if r == "FAIL":
            bad.append(l.strip())
        else:
            classes[tag.split(".")[0][2:]][l.split(": ")[0].split(" ", 2)[2]] += 1
    L = leaves.get(tag)
    if end is None or end != L or n != len(range(c, L or 0, k)):
        bad.append(f"EA chunk {ch}: {n} records, END {end}, expected {len(range(c, L or 0, k))} of {L} leaves")
for tag, h in hashes.items():
    if any(v > 1 for v in h.values()):
        bad.append(f"{tag}: duplicate leaf hashes in the EA records")
if sum(res.values()) != sum(leaves.values()):
    bad.append(f"EA records {sum(res.values())} != leaves of the search jobs {sum(leaves.values())}")
print(f"search jobs: {agree} agree, {differ} differ (of 67); EA chunks done {len(have_done)}, ended {ends}; leaves equivalent "
      f"{res['equivalent']}, failed {res['FAIL']}")
for t in sorted(classes):
    print(f"  type {t}: " + "; ".join(f"{c} {n}" for c, n in sorted(classes[t].items())))
for b in bad[:20]:
    print("  !", b)
ok = not bad and agree == 67 and differ == 0 and len(have_done) == len(want_chunks) == 641
print("FPGAP CHECK OK" if ok else "FPGAP CHECK FAILED")
sys.exit(0 if ok else 1)
