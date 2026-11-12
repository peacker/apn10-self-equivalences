#!/usr/bin/env python3
"""Check (or write) MANIFEST.sha256 at the repository root: SHA-256 of every tracked file except the manifest.

    python3 -I tools/check_manifest.py           # check every listed file; exit status 1 on any mismatch
    python3 -I tools/check_manifest.py --write   # rewrite the manifest from `git ls-files` (maintainers)

The format is that of `shasum -a 256` / `sha256sum`, so `shasum -a 256 -c MANIFEST.sha256` gives the same check."""
import hashlib, os, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAN = os.path.join(ROOT, "MANIFEST.sha256")


def sha(path):
    h = hashlib.sha256()
    with open(os.path.join(ROOT, path), "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if "--write" in sys.argv:
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, check=True).stdout
    files = sorted(f for f in files.decode().split("\0") if f and f != "MANIFEST.sha256")
    with open(MAN, "w") as out:
        for f in files:
            out.write("%s  %s\n" % (sha(f), f))
    print("wrote %d entries" % len(files))
    sys.exit(0)

bad = 0; n = 0
for line in open(MAN):
    digest, path = line.rstrip("\n").split("  ", 1)
    n += 1
    if not os.path.isfile(os.path.join(ROOT, path)):
        print("MISSING", path); bad += 1
    elif sha(path) != digest:
        print("CHANGED", path); bad += 1
print("%d files checked, %d mismatches" % (n, bad))
sys.exit(1 if bad else 0)
