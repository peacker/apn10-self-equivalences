#!/usr/bin/env python3
"""Type 47: check from the stored job outputs that every (class, X) slice of stage 2 was searched to the end
and that no slice gave a function.

1. stage 1 / rank-h lemma: `stage2i` with LISTP=1 on each of the 806 classes in cls/ reports "dead on U x U"
   (the rank-h lemma kills the class) or the list of X representatives grouped by prefix. Expected: 703 dead,
   103 alive, 8,056 (class, X) slices in total, exactly the slices of jobs47n.txt.
2. finished jobs: results/t47_out.tar.gz (unpacked to results/out/) holds the outputs of the stage2g runs
   (TAG = CLASS.X) and of the stage2i runs (TAG = CLASS.X or CLASS.pG, G a prefix group). A job counts if
   TAG.done says "rc 0" and TAG.out has the summary line "... found 0 (APN 0)".
3. coverage: the union of the slices of the finished jobs must be all 8,056 slices (and equal
   results/covered_u.txt). Also prints the summed run time of the counted jobs.
Job records that appear more than once (the stage2g outputs were synced from two servers) must agree: same .out and
same return code, else the check fails. The LISTP runs must exit with status 0. Exit status 1 on any failure.
usage: check_coverage.py [--stage2i PATH]   (default: compile stage2i.c afresh into a temporary directory)"""
import collections, glob, os, re, shutil, subprocess, sys, tarfile, tempfile
here = os.path.dirname(os.path.abspath(__file__))
res = os.path.join(here, "results")
out = os.path.join(res, "out")
if not os.path.isdir(out):
    with tarfile.open(os.path.join(res, "t47_out.tar.gz")) as t:
        try:
            t.extractall(out, filter="data")
        except TypeError:
            t.extractall(out)
tmp = tempfile.mkdtemp(prefix="t47cov_")      # inside the caller's TMPDIR, if set
# no solver options inherited; TMPDIR points to our own directory so the compiler's temporary files stay there
ENV = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "LC_ALL": "C", "TMPDIR": tmp}
if "--stage2i" in sys.argv:
    exe = sys.argv[sys.argv.index("--stage2i") + 1]
else:
    exe = os.path.join(tmp, "stage2i")
    subprocess.run(["cc", "-O3", "-DKU=3", "-DMZ=4", "-o", exe, os.path.join(here, "stage2i.c")], check=True, env=ENV)
errors = []

groups, dead = {}, []
for f in sorted(glob.glob(os.path.join(here, "cls", "*.txt"))):
    c = os.path.basename(f)[:-4]
    p = subprocess.run([exe, "cls/" + c + ".txt"], cwd=here, env=dict(ENV, LISTP="1"), capture_output=True, text=True)
    r = p.stdout
    if p.returncode != 0:
        errors.append(f"LISTP run on {c} exited with status {p.returncode}"); continue
    if "dead on U x U" in r:
        dead.append(c); continue
    g = {}
    for line in r.splitlines():
        m = re.match(r"group (\d+):(.*)", line)
        if m: g[int(m.group(1))] = [int(x) for x in m.group(2).split()]
    if not g or not re.search(r"groups \d+, X \d+", r):
        errors.append(f"LISTP run on {c}: no group list")
    groups[c] = g
if tmp:
    shutil.rmtree(tmp, ignore_errors=True)
slices = {(c, x) for c, g in groups.items() for xs in g.values() for x in xs}
jobs = {(l.split()[0], int(l.split()[1])) for l in open(os.path.join(here, "jobs47n.txt"))}
print(f"classes: {len(dead) + len(groups)}, dead on U x U {len(dead)}, alive {len(groups)}, slices {len(slices)}"
      f" (jobs47n.txt: {len(jobs)}, same set: {slices == jobs})")

copies = collections.defaultdict(list)
for d in sorted(glob.glob(os.path.join(out, "**", "*.done"), recursive=True)):
    copies[os.path.basename(d)[:-5]].append(d)
covered, secs, njobs, bad, ndup, disagree = set(), 0, 0, [], 0, []
for tag, ds in sorted(copies.items()):
    recs = []
    for d in ds:
        o = d[:-5] + ".out"
        recs.append((open(d).read().split()[:2], open(o).read() if os.path.exists(o) else None))
    if len(ds) > 1:
        ndup += 1
        if any(x != recs[0] for x in recs[1:]):
            disagree.append(tag)
    rc, summ = recs[0][0], recs[0][1] or ""
    c, rest = tag.rsplit(".", 1)
    if rc != ["rc", "0"] or "found 0 (APN 0)" not in summ or not summ.startswith(f"cls/{c}.txt:") \
            or "limit" in summ.lower():
        bad.append(tag); continue
    if c not in groups:
        errors.append(f"record {tag}: class not alive"); continue
    if rest.startswith("p"):
        g = int(rest[1:])
        if g not in groups[c] or f"prefix {g} ({len(groups[c][g])} X)" not in summ:
            errors.append(f"record {tag}: prefix group does not match the LISTP output"); continue
        xs = groups[c][g]
    else:
        xs = [int(rest)]
        if (c, xs[0]) not in slices:
            errors.append(f"record {tag}: X not a representative"); continue
    covered |= {(c, x) for x in xs}
    secs += int(re.search(r"secs (\d+)", open(ds[0]).read()).group(1)); njobs += 1
if disagree:
    errors.append(f"{len(disagree)} duplicated tags disagree between copies, e.g. {disagree[:3]}")
cu = {(l.strip().rsplit(".", 1)[0], int(l.strip().rsplit(".", 1)[1])) for l in open(os.path.join(res, "covered_u.txt"))}
print(f"finished jobs counted: {njobs} distinct tags (rejected {len(bad)}; {ndup} tags stored more than once, all copies "
      f"agree: {not disagree}), slices covered {len(covered & slices)} of {len(slices)},"
      f" equal to covered_u.txt: {covered == cu}; summed job time {secs / 3600:.0f} h (wall time per job, all jobs found 0)")
for e in errors[:20]:
    print("  !", e)
ok = (not errors and len(dead) == 703 and len(groups) == 103 and slices == jobs and len(slices) == 8056
      and slices <= covered and covered == cu)
print("TYPE 47 COVERAGE OK" if ok else "TYPE 47 COVERAGE FAILED")
sys.exit(0 if ok else 1)
