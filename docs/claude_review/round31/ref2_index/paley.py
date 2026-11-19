import sys
from itertools import combinations
from lat import lattice_index, anchored_index
def chi(x,q):
    x%=q
    if x==0: return 0
    return 1 if pow(x,(q-1)//2,q)==1 else -1
def core(q):
    M=q+1
    # marks: 0 = border, 1..q = field elements 0..q-1
    K=[[0]*M for _ in range(M)]
    for x in range(M):
        for y in range(M):
            if x==0 and y==0: K[x][y]=1
            elif x==0: K[x][y]=-1
            elif y==0: K[x][y]=1
            else: K[x][y]=(1 if x==y else 0)+chi((y-1)-(x-1),q)
    cols=[]
    for a in range(1,M):
        s1=[K[x][a] for x in range(M)]
        s2=[K[a][x]*K[0][x] for x in range(M)]  # K^T(x,a)K^T(x,0)
        cols.append(s1); cols.append(s2)
    # check K K^T = M I
    for x in range(M):
        for y in range(M):
            assert sum(K[x][t]*K[y][t] for t in range(M))==(M if x==y else 0)
    B=[[(1+s[x])//2 for s in cols] for x in range(M)]  # rows = marks
    return B
for q in [3,7,11,19,23,31]:
    M=q+1
    B=core(q)
    # check row separation in Hamming (all primes equal weight): pair distance
    W=2*q
    dmin=min(sum(abs(B[i][j]-B[k][j]) for j in range(W)) for i,k in combinations(range(M),2))
    full=anchored_index(B,tuple(range(M)))
    bad5=sum(1 for T in combinations(range(M),5) if anchored_index(B,T)!=1)
    n5=sum(1 for _ in combinations(range(M),5))
    line=f"q={q} M={M} primes={W} min pair Hamming={dmin} (W/2={W/2}) whole index={full} (M/2={M//2}) five-subsets with index!=1: {bad5}/{n5}"
    if M<=12:
        badproper=0; tot=0
        for k in range(2,M):
            for T in combinations(range(M),k):
                tot+=1
                if anchored_index(B,T)!=1: badproper+=1
        line+=f"; proper subsets (size>=2) with index!=1: {badproper}/{tot}"
    print(line, flush=True)
