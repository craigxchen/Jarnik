"""Cut profile of Astra's mixed-rate four-point Pell example (rates 6,7,7) and equal-rate (1,1,1)."""
import math, itertools
from glib import *
from sympy import factorint
def adv(v,t):
    x,y=v
    for _ in range(t): x,y=5*x+2*y,2*x+y
    return (x,y)
def run(rates,n):
    seeds=((3,1),(4,1),(5,1)); G=[adv(v,r*n) for v,r in zip(seeds,rates)]
    pts=[]
    for fl in ((),(0,1),(0,2),(1,2)):
        p=(1,0) if not fl else (0,-1)
        for j,z in enumerate(G): p=mul(p,z if j in fl else conj(z))
        pts.append(p)
    g=pts[0]
    for p in pts[1:]: g=gcd(g,p)
    pts=[exdiv(p,g) for p in pts]; H=nrm(pts[0]); W=math.log(H)
    fac={}
    for z in G:
        for p,e in factorint(nrm(z)).items(): fac[p]=1
    cuts={}
    for p in fac:
        if p%4!=1: continue
        pi=split_prime(p); a=[val(x,pi) for x in pts]; e=max(a)+0
        b=[val(x,conj(pi)) for x in pts]; e=a[0]+b[0]
        for t in range(1,e+1):
            hi=frozenset(k for k in range(4) if a[k]>=t)
            if len(hi) in (0,4): continue
            side=hi if len(hi)<2 or (len(hi)==2 and 0 in hi) else frozenset(range(4))-hi
            cuts[side]=cuts.get(side,0)+math.log(p)
    d2=max(nrm((x[0]-y[0],x[1]-y[1])) for x,y in itertools.combinations(pts,2)); ld=(math.log(d2)-W)/2
    print('rates',rates,'n',n,'logH=%.1f'%W,'width exponent %.4f'%(-ld/W),
          {tuple(sorted(c)):round(w/W,3) for c,w in cuts.items()}, 'Wtot/W=%.3f'%(sum(cuts.values())/W))
for n in (1,2): run((6,7,7),n)
for n in (3,6): run((1,1,1),n)
