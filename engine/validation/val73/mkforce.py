#!/usr/bin/env python3
"""-force file with beta(e_i, e_j), i < j < U, of the LUT on line LINE of LUTFILE.  usage: mkforce.py LUTFILE LINE U"""
import sys
F = [list(map(int, l.split(","))) for l in open(sys.argv[1]) if "," in l][int(sys.argv[2])]
U = int(sys.argv[3])
for j in range(U):
    for i in range(j):
        print(j * (j - 1) // 2 + i, F[(1 << i) | (1 << j)] ^ F[1 << i] ^ F[1 << j] ^ F[0])
