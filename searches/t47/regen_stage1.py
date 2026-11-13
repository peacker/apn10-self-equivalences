#!/usr/bin/env python3
"""Type 47 stage 1: regenerate the 806 (S, V) class representatives with stage1.py and stage1b.py in a temporary
directory and compare them with the stored ones (cls/*.txt byte for byte, cls_summary.txt, the contents of aniso_classes.pkl).

stage1.py enumerates the anisotropic Hermitian spans S (dim 3 and 4) and their orbits under GL_3(F_4) (union-find
with 4 random generators); stage1b.py computes, for each orbit representative, Stab(S) by brute force over all
181,440 elements of GL_3(F_4) and asserts |Stab(S)| * |orbit| = |GL_3(F_4)| (so the orbits are complete), then the
orbits of Stab(S) on the F_4-subspaces V and the -force file of each (S, V). Both are deterministic (fixed seed).
About 40 s. usage: python3 -I searches/t47/regen_stage1.py   (exit status 1 on any difference)"""
import filecmp, glob, os, pickle, shutil, subprocess, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp(prefix="t47stage1_")
env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "LC_ALL": "C", "TMPDIR": tmp}
ok = True
try:
    for s in ("stage1.py", "stage1b.py"):
        r = subprocess.run([sys.executable, "-I", os.path.join(here, s)], cwd=tmp, env=env, capture_output=True, text=True)
        if s == "stage1b.py":
            open(os.path.join(tmp, "cls_summary.txt"), "w").write(r.stdout)
        print(f"{s}: exit {r.returncode}; " + " | ".join(r.stderr.strip().splitlines()[-2:]))
        ok &= r.returncode == 0
    new = sorted(os.path.basename(f) for f in glob.glob(os.path.join(tmp, "cls", "*.txt")))
    old = sorted(os.path.basename(f) for f in glob.glob(os.path.join(here, "cls", "*.txt")))
    same_names = new == old
    _, mismatch, errs = filecmp.cmpfiles(os.path.join(tmp, "cls"), os.path.join(here, "cls"), old, shallow=False)
    same_sum = filecmp.cmp(os.path.join(tmp, "cls_summary.txt"), os.path.join(here, "cls_summary.txt"), shallow=False)
    load = lambda d: pickle.load(open(os.path.join(d, "aniso_classes.pkl"), "rb"))     # compare contents, not bytes
    same_pkl = load(tmp) == load(here)
    print(f"classes regenerated {len(new)}, stored {len(old)}; same names {same_names}; differing files "
          f"{len(mismatch) + len(errs)}; cls_summary.txt equal {same_sum}; aniso_classes.pkl equal {same_pkl}")
    ok &= same_names and not mismatch and not errs and same_sum and same_pkl and len(new) == 806
finally:
    shutil.rmtree(tmp, ignore_errors=True)
print("TYPE 47 STAGE 1 REPRODUCED" if ok else "TYPE 47 STAGE 1 DIFFERS")
sys.exit(0 if ok else 1)
