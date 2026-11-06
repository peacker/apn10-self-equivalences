# Monte-Carlo: fraction of random 6-dim F_2-subspaces of Herm_4(F_4) with every nonzero element of rank >= 3
import random
F4M=[[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]; CONJ=[0,1,3,2]; INV=[0,1,3,2]
def mat(v):
    H=[[0]*4 for _ in range(4)]; k=0
    for i in range(4): H[i][i]=(v>>k)&1; k+=1
    for i in range(4):
        for j in range(i+1,4):
            h=(v>>k)&3; k+=2; H[i][j]=h; H[j][i]=CONJ[h]
    return H
def rank(M):
    M=[r[:] for r in M]; r=0
    for c in range(4):
        p=next((i for i in range(r,4) if M[i][c]),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]; iv=INV[M[r][c]]; M[r]=[F4M[iv][x] for x in M[r]]
        for i in range(4):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[M[i][k]^F4M[f][M[r][k]] for k in range(4)]
        r+=1
    return r
RK=[rank(mat(v)) for v in range(1<<16)]
cnt={}
for v in range(1,1<<16): cnt[RK[v]]=cnt.get(RK[v],0)+1
print("rank counts in Herm_4:",cnt)
random.seed(1); good=0; T=200000
for _ in range(T):
    b=[]
    while len(b)<6:
        v=random.getrandbits(16)
        x=v
        for y in b: x=min(x,x^y)
        if x: b.append(x)
    ok=True
    for m in range(1,64):
        v=0
        for k in range(6):
            if (m>>k)&1: v^=b[k]
        if RK[v]<3: ok=False; break
    good+=ok
p=good/T
# number of 6-dim subspaces of F_2^16
num=1;den=1
for i in range(6): num*=(2**16-2**i); den*=(2**6-2**i)
G=255*252*240*192//3*2
print("fraction %.4g, subspaces %.3g, good %.3g, group %.3g, orbits >= %.3g"%(p,num/den,p*num/den,G,p*num/den/G))
