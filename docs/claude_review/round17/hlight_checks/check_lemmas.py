# Exact checks of the lemmas used in hlight.md (all integer arithmetic).
# L1 identities; L2 zero bound Z(f*chi) <= (q-1)/2 - j0 for 1<=|supp f|<=(q+1)/2, ||f||_inf<q;
# L3 case bounds: |F|=1 ternary; F=0 ternary sigma >= max(n-1, g(n)) (n>=4), g(n)=(q(n-2)+4n-n^2)/(n+2).
import numpy as np, sys
from fractions import Fraction as Fr
from hlib import *
rng=np.random.default_rng(7)
def j0(f,q):
    T=np.nonzero(f)[0]
    for j in range(0,(q-1)//2+1):
        if sum(int(f[t])*pow(int(t),j,q) for t in T)%q: return j
    return None
def g(n,q): return Fr(q*(n-2)+4*n-n*n,n+2)
viol={'Z':0,'F1':0,'F0':0,'id':0}; tight={'Z':0}
for q in [7,11,19,23,31,43,47,59,67,71,79,83,103,127,199]:
    ch=chi_vec(q); cnt=0
    for trial in range(1500 if q<100 else 400):
        # random ternary or small-integer f with prescribed support size
        n=int(rng.integers(1,(q+1)//2+1))
        f=np.zeros(q,dtype=np.int64); T=rng.choice(q,size=n,replace=False)
        if trial%3==0: f[T]=rng.choice([-2,-1,1,2],size=n)
        else: f[T]=rng.choice([-1,1],size=n)
        if trial%2==0 and n>=2:  # force F=0 for ternary
            if np.all(np.abs(f)<=1):
                k=n//2; f[T[:k]]=1; f[T[k:2*k]]=-1; f[T[2*k:]]=0
        if not f.any(): continue
        c=conv_chi(f,q,ch); Z=int((c==0).sum()); J=j0(f,q)
        if J is None or Z>(q-1)//2-J: viol['Z']+=1
        if Z==(q-1)//2-J: tight['Z']+=1
        F=int(f.sum()); lam=np.append(f,-F); s=sigma_formula(f,q,ch)
        # identity sum t^2
        U=(F+f+c)//2; V=(F+f-c)//2
        if int((U*U+V*V).sum())*2!=(q+1)*int((lam*lam).sum()): viol['id']+=1
        nn=int(np.abs(lam).sum())
        tern=np.all(np.abs(f)<=1)
        if tern and abs(F)==1 and nn>=4:
            bound=Fr(q-1,2) if nn==4 else Fr(q+1,2)
            if s<bound: viol['F1']+=1
        if tern and F==0 and nn>=4:
            bound=max(Fr(nn-1),g(nn,q)) if nn<=(q+1)//2 else Fr(nn-1)
            if s<bound: viol['F0']+=1
        cnt+=1
    print(f"q={q} tested={cnt} cumulative violations={viol} tight-zero-bound cases={tight['Z']}",flush=True)
# pairs attain zero bound
for q in [19,43,199]:
    ch=chi_vec(q); f=np.zeros(q,dtype=np.int64); f[0]=1; f[5]=-1
    c=conv_chi(f,q,ch); print(f"q={q} pair: zeros={(c==0).sum()} bound (q-3)/2={(q-3)//2}")
