#!/usr/bin/env python3
"""Master table of the 299 prime-order self-equivalence types at n = 10: how each is settled (proof, search, both),
search cost, and the quadratic APN classes it contains.
Inputs: types10.tsv, proofs10.txt (criteria C/R/Q/W=0), lemma_dual.txt, bbf_excluded.txt (Baudrin et al. 2025 bound), search_status.txt (per-type search
completeness from the server .done files, 2026-11-05), server_done.txt (first line of every .done file, with secs;
machine labels srv1, srv2 replace the original server names).
Special runs (t24, t47, t73, t41s, type 64) are entered by hand below. Output: master_table.tsv (or OUT).
usage: master_table.py [-o OUT]"""
import collections, os, re, sys
here = os.path.dirname(os.path.abspath(__file__))
P = lambda f: os.path.join(here, f)
types = [l.rstrip("\n").split("\t") for l in open(os.path.join(here, "..", "types", "types10.tsv"))]
proof = {}
for l in open(P("proofs10.txt")):
    t = l.rstrip("\n").split("\t")
    c = t[6]
    if c != "-": proof[int(t[0])] = " ".join(x.split("(")[0] for x in c.split()) if c != "W=0" else "W=0"
for t in (130, 169): proof[t] = "Wong"
for t in (55, 56): proof[t] = "lemma56"
dual = {int(l.split()[1]) for l in open(P("lemma_dual.txt")) if l.startswith("type ")}
for t in dual:
    proof[t] = (proof[t] + " D") if t in proof else "D"
bbf = {int(x) for x in open(P("bbf_excluded.txt")).read().split()}
for t in bbf:
    proof[t] = (proof[t] + " BBF") if t in proof else "BBF"
status = {}
for l in open(P("search_status.txt")):
    t, rest = l.split(" ", 1)
    status[int(t)] = rest.strip()
secs = collections.defaultdict(float)                      # (machine, dir, type) -> seconds
for l in open(P("server_done.txt")):
    m = re.match(r"(\S+) (\S+)/t_(\d+)(?:\.\d+)?\.done .*secs (\d+)", l)
    if m: secs[(m.group(1), m.group(2), int(m.group(3)))] += int(m.group(4))
def search_of(t):
    s = status.get(t)
    if not s: return "none", "", ""
    best = None
    for part in s.split(" | "):
        port_dir, st, sh, lv = part.split()
        port, d = port_dir.split(":")
        if st == "COMPLETE": best = (port, d, lv.split("=")[1])
    if not best: return "partial", "", ""
    return "complete", "%.1f" % (secs[(best[0], best[1], t)] / 3600), best[2]
CLASSES = {
    24: "Gold x^3; Gold x^9; F4 a=1",
    73: "24 classes: 20 known + 4 new (see PAPER_NOTES)",
    64: "16 known (conjectured complete)",
    129: "Gold x^3; Gold x^9; Beierle-Leander #5; F13 k=2", 170: "F3 i=1 c=1014", 171: "F3 i=2 c=941", 175: "Gold x^3", 180: "Gold x^9",
    176: "Beierle-Leander #3; F12 x6; F13 x4", 231: "Beierle-Leander #1"}
SPECIAL = {  # type -> (search, core-hours, leaves)
    24: ("complete (t24)", "7500", "3"), 47: ("complete (t47)", "13261", "0"),
    73: ("complete (sym10 + rerun)", "12600", "4536"), 41: ("partial (t41s; empty by BBF, C, D)", "", "0"),
    64: ("not exhaustive (4 searches)", "", "")}
rows = []; summary = collections.Counter()
for t in types:
    ty = int(t[0]); p = t[1]; kind = t[2]
    s, ch, lv = SPECIAL.get(ty) or search_of(ty)
    pr = proof.get(ty, "-")
    cl = CLASSES.get(ty, "-")
    if ty == 64: how = "open"
    elif pr != "-" and s.startswith("complete"): how = "both"
    elif pr != "-": how = "proof"
    elif s.startswith("complete"): how = "search"
    else: how = "UNSETTLED"
    summary[how] += 1
    if cl != "-": summary["contains APN"] += 1
    rows.append((ty, p, kind, t[3], t[4], how, pr, s, ch, lv, cl))
outp = sys.argv[sys.argv.index("-o") + 1] if "-o" in sys.argv else P("master_table.tsv")
with open(outp, "w") as f:
    f.write("type\tp\tkind\tA\tB\tsettled_by\tproof\tsearch\tcore_hours\tleaves\tclasses\n")
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
print(dict(summary))
for r in rows:
    if r[5] in ("UNSETTLED", "open") or r[10] != "-": print("\t".join(map(str, r)))
