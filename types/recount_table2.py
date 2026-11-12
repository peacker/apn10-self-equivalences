#!/usr/bin/env python3
# Independent recount of Table II of the paper (299 types). Self-contained, does not use types/enumerate.py.
# Output: one line per p: p, #(p,p), #(p,1), #(1,p), total; then "total 299".
# usage: python3 -I types/recount_table2.py
# Independent recount of Table II: prime-order types in GL(10,2) up to independent conjugacy and simultaneous powering.
from itertools import product
n=10
def ordp(p):
    k,x=1,2%p
    while x!=1: x=x*2%p; k+=1
    return k
def factors(p):
    # irreducible factors of Phi_p over F2 <-> cyclotomic cosets of 2 mod p (excluding 0)
    seen=set(); cos=[]
    for a in range(1,p):
        if a in seen: continue
        c=set(); x=a
        while x not in c: c.add(x); x=x*2%p
        seen|=c; cos.append(frozenset(c))
    return cos
total=0
for p in [2,3,5,7,11,17,31,73,127]:
    if p==2:
        classes=list(range(1,6)); act=lambda r,c:c; rs=[1]
    else:
        cos=factors(p); d=ordp(p); m=len(cos)
        # class = tuple of multiplicities, total degree <= n, not all zero
        classes=[t for t in product(range(n//d+1),repeat=m) if 0<sum(t)*d<=n]
        def act(r,t,cos=cos):
            # power r maps factor with root z to root z^r: coset C -> r*C
            idx={c:i for i,c in enumerate(cos)}
            out=[0]*len(cos)
            for i,c in enumerate(cos): out[idx[frozenset(r*x%p for x in c)]]+=t[i]
            return tuple(out)
        rs=list(range(1,p))
    def orbits(items,f):
        seen=set(); k=0
        for it in items:
            if it in seen: continue
            k+=1
            for r in rs: seen.add(f(r,it))
        return k
    one=orbits(classes,lambda r,c:act(r,c))
    pp=orbits([(a,b) for a in classes for b in classes],lambda r,ab:(act(r,ab[0]),act(r,ab[1])))
    print(p,pp,one,one,pp+2*one); total+=pp+2*one
print("total",total)
