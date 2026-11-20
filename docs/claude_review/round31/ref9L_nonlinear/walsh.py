# Check (IX.5) coefficient bound and integrality for bits that are functions of <= d affine forms on F_2^t.
from itertools import product, combinations
from math import isqrt
import sys
def dot(u,x): return bin(u&x).count('1')&1
def juntas(t,d):
    # all Boolean functions on F_2^t that factor through <= d affine forms
    seen=set()
    forms=list(range(1<<t))  # linear parts; affine constant irrelevant up to relabel of g
    for ls in combinations(forms,d) if d>0 else [()]:
        for g in range(1<<(1<<d)):
            for c in range(1<<d):
                tab=tuple((g>>(sum(((dot(l,x)^((c>>k)&1))<<k) for k,l in enumerate(ls))))&1 for x in range(1<<t))
                seen.add(tab)
    return seen
def rank2(vs):
    rows=[];
    for v in vs:
        for r in rows: v=min(v,v^r)
        if v: rows.append(v)
    return len(rows)
def check(t,d):
    m0=2**(t-d); fs=juntas(t,d)
    nz=list(range(1,1<<t))
    best=0; nonint=0
    # all bases U of F_2^t
    bases=[U for U in combinations(nz,t) if rank2(list(U))==t]
    for f in fs:
        S={u:sum((-1)**dot(u,x)*f[x] for x in range(1<<t)) for u in nz}
        for u in nz:
            if S[u]%m0: nonint+=1
        # max over bases of sum |S/m0|
        for U in bases:
            v=sum(abs(S[u]) for u in U)//m0
            if v>best: best=v
    print(f"t={t} d={d} m0={m0} #functions={len(fs)} #bases={len(bases)} nonintegral={nonint} max sum_U|b|={best}")
check(2,1); check(4,2); check(3,2)
# d=3 via quotient formula: b = 2^{d-d'} S on F_2^{d'}; max over d'<=3, g, independent u's (<=d' of them in K^perp)
for d in (2,3,4):
    best=0
    for dp in range(0,d+1):
        nz=list(range(1,1<<dp))
        for g in range(1<<(1<<dp)):
            S=[sum((-1)**dot(u,x)*((g>>x)&1) for x in range(1<<dp)) for u in range(1<<dp)]
            for k in range(0,dp+1):
                for U in combinations(nz,k):
                    if rank2(list(U))==k:
                        best=max(best,2**(d-dp)*sum(abs(S[u]) for u in U))
    print("d",d,"max sum",best,"claimed b_d",int((2**(d-1))*d**0.5) if d!=4 else 16)
