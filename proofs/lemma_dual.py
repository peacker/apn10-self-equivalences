#!/usr/bin/env python3
"""Dual intertwiner lemma (2026-11-05). Let beta(Ax, Ay) = B beta(x, y), A of odd prime order, F APN. For z in Fix(A)
and u in the A-isotypic part U_phi = ker phi(A) (phi an irreducible factor of x^p - 1, phi != x + 1), beta(Au, z) =
beta(Au, Az) = B beta(u, z), so u -> beta(u, z) intertwines A and B and beta(U_phi, Fix A) lies in ker phi(B).
For u in U_phi - {0} the map z -> beta(u, z) on Fix(A) has kernel in {0, u} ∩ Fix(A) = {0}, so it is injective:
dim Fix(A) <= dim ker phi(B) whenever U_phi != 0. Prints the n = 10 types this excludes (types10.tsv).
usage: lemma_dual.py"""
import os, sys
n=10
def app(M,v):
    r=0
    for i in range(n):
        if (v>>i)&1: r^=M[i]
    return r
def polyeval(coefs,M):  # coefs: list of bits, c0 + c1 x + ...; returns matrix columns
    cols=[0]*n
    for i in range(n):
        v=1<<i; acc=0; p=v
        for c in coefs:
            if c: acc^=p
            p=app(M,p)
        cols[i]=acc
    return cols
def kerdim(cols):
    rows=list(cols); r=0  # rank of the map = rank of column set
    basis=[]
    for v in rows:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    return n-len(basis)
# irreducible factors of x^p - 1 over F2 for small p, as bit lists (low degree first)
def polymulmod2(a,b):
    r=0
    while b:
        if b&1: r^=a
        a<<=1; b>>=1
    return r
def polydivmod(a,b):
    q=0; db=b.bit_length()
    while a.bit_length()>=db:
        s=a.bit_length()-db; q|=1<<s; a^=b<<s
    return q,a
def factors(p):
    f=(1<<p)|1; out=[]
    d=2
    while f.bit_length()>1:
        found=False
        for g in range(1<<d, 1<<(d+1)):
            if polydivmod(f,g)[1]==0 and g.bit_length()-1==d:
                # check irreducible by trial
                irr=all(polydivmod(g,h)[1]!=0 for h in range(2,g) if 0<h.bit_length()-1<d)
                if irr:
                    out.append(g); f=polydivmod(f,g)[0]; found=True; break
        if not found: d+=1
        if d>12: break
    return out
rows=[l.rstrip("\n").split("\t") for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "types", "types10.tsv"))]
excl=[]
for r in rows:
    tid,p=r[0],int(r[1]); A=list(map(int,r[5].split(","))); B=list(map(int,r[6].split(",")))
    if p==2: continue
    fixA=kerdim([app(A,1<<i)^(1<<i) for i in range(n)])
    if fixA==0: continue
    for g in factors(p):
        if g==3: continue  # x+1
        bits=[(g>>i)&1 for i in range(g.bit_length())]
        dA=kerdim(polyeval(bits,A)); dB=kerdim(polyeval(bits,B))
        if dA>0 and fixA>dB:
            excl.append((tid,p,fixA,bin(g),dA,dB)); break
for e in excl: print("type %s p=%d dimFixA=%d phi=%s dimA_phi=%d dimB_phi=%d"%e)
print(len(excl),"types excluded by the dual check")
