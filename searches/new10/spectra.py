#!/usr/bin/env python3
"""Ortho-derivative differential and extended Walsh spectra of the n = 10 candidates (one LUT per line)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
from apnsearch.invariants import ortho_invariant  # noqa: E402
for line in open(sys.argv[1]):
    F = np.array(list(map(int, line.split(","))), dtype=np.int64)
    d, w = ortho_invariant(F, 10)
    print("diff:", dict(d)); print("walsh:", dict(w)); print()
