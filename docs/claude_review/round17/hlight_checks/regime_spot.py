# Spot check at regime q=3251: Lemma 1.1 formula vs direct pattern sums, and sigma >= LB(q)=1527.5 on random non-pair lambda.
import numpy as np
from hlib import *
q=3251; P,Pp=patterns(q); ch=chi_vec(q); rng=np.random.default_rng(3); bad=0; mism=0; mn=10**9
P=[np.array(S) for S in P]; Pp=[np.array(S) for S in Pp]
for trial in range(60):
    n=int(rng.choice([3,4,5,6,8,20,200,1600]))
    f=np.zeros(q,dtype=np.int64); T=rng.choice(q,size=n,replace=False); f[T]=rng.choice([-1,1],size=n)
    if trial%2==0: k=n//2; f[T[:k]]=1; f[T[k:2*k]]=-1; f[T[2*k:]]=0
    lam=np.append(f,-f.sum())
    if np.abs(lam).sum()<=2 or np.abs(lam).max()==0: continue
    s1=sigma_formula(f,q,ch)
    if trial<12:
        s2=sigma_direct(lam,q,P,Pp); mism+= (s1!=s2)
    mn=min(mn,s1); bad+= (s1<1527.5)
print(f"q={q}: formula-vs-direct mismatches={mism} (12 checked), violations of LB={bad}, min sigma seen={mn}")
