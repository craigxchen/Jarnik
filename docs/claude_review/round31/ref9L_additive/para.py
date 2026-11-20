#!/usr/bin/env python3
"""Minimal normalised angular diameter of 'signature parallelograms' on genuine circles:
4 distinct points in one unit class whose unary signatures satisfy sig_i^sig_j = sig_k^sig_l
(pairs {i,j},{k,l} disjoint). IX.3 circle lemma (recomputed constant) predicts
diam * N'^{1/4} >= 2 where N' is the normalised norm of the 4 points."""
import math, itertools, sys
from collections import defaultdict
from circle import points
from sympy import primerange
split=[p for p in primerange(5,200) if p%4==1]
facs=[{p:1 for p in split[:k]} for k in range(3,10)]+[{p:2 for p in split[:k]} for k in range(2,6)]+[{5:3,13:3,17:2,29:1},{5:5,13:2,17:1,29:1,37:1}]
best=(1e9,None)
for fac in facs:
    ps,pts=points(fac)
    for u in range(4):
        P=[(t,r) for t,uu,r,w in pts if uu==u]
        def usig(r):
            v=0;off=0
            for k,p in enumerate(ps):
                v|=((1<<r[k])-1)<<off; off+=fac[p]
            return v
        S=[usig(r) for t,r in P]
        g=defaultdict(list)
        for i,j in itertools.combinations(range(len(P)),2): g[S[i]^S[j]].append((i,j))
        loc=1e9
        for h,L in g.items():
            if len(L)<2: continue
            # for each pair of disjoint pairs compute normalised diameter; limit work
            for (i,j),(k,l) in itertools.combinations(L,2):
                Q=[i,j,k,l]
                if len(set(Q))<4: continue
                ang=sorted(P[x][0] for x in Q)
                gaps=[ang[(a+1)%4]-ang[a] for a in range(4)]; gaps[3]+=2*math.pi
                diam=2*math.pi-max(gaps)
                Np=1
                for kk,p in enumerate(ps):
                    vals=[P[x][1][kk] for x in Q]; Np*=p**(max(vals)-min(vals))
                val=diam*Np**0.25
                if val<loc: loc=val; ex=(fac,Q,[P[x][1] for x in Q],Np)
        if loc<best[0]: best=(loc,ex)
    print(fac, 'min diam*N\'^(1/4) so far %.4f'%best[0], flush=True)
print('BEST',best)
