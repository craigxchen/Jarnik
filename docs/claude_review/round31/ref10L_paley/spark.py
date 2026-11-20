# spark of T=[1; S_M] for M = complement of half-row set {i: K_i0 = 1}; full column rank of all (k+1)-column subsets?
import sys, random
from itertools import combinations
from paley1 import chi, S_of
PR=2**61-1
def rank_mod(A):
    A=[[x%PR for x in r] for r in A]; rk=0; rows=len(A); cols=len(A[0])
    for c in range(cols):
        p=next((i for i in range(rk,rows) if A[i][c]),None)
        if p is None: continue
        A[rk],A[p]=A[p],A[rk]; inv=pow(A[rk][c],-1,PR)
        for i in range(rows):
            if i!=rk and A[i][c]:
                f=A[i][c]*inv%PR; A[i]=[(A[i][t]-f*A[rk][t])%PR for t in range(cols)]
        rk+=1
    return rk
def test(m,exhaustive,samples=20000):
    S=S_of(m); half=[i for i in range(m) if chi(i-0,m)==1]  # K_i0 = chi(0-i)
    half=[i for i in range(m) if chi(0-i,m)==1]
    M=[i for i in range(m) if i not in half]; k=len(M)
    T=[[1]*m]+[S[i] for i in M]
    it=combinations(range(m),k+1) if exhaustive else (random.sample(range(m),k+1) for _ in range(samples))
    bad=0;cnt=0
    for U in it:
        cnt+=1
        if rank_mod([[T[r][c] for c in U] for r in range(k+1)])<k+1: bad+=1
    print("m=%d retained=%d k=%d: (k+1)-column subsets tested %d, rank-deficient (mod p) %d"%(m,len(half),k,cnt,bad))
random.seed(2)
test(11,True); test(19,True); test(23,False); test(43,False,3000)
