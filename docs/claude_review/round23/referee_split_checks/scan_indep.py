"""Independent low-height cross-check of split.md scan.c: all primitive 5-tuples of distinct positive integers
with lcm = N <= NMAX and gamma = log H/(-log dist) <= 5, H=(max/min)*N, dist=(max-min)/min.
Different enumeration: sieve divisor lists, then sliding window of ratio <= 1+N^(-1/5)."""
import sys, math, itertools
from math import gcd
NMAX=int(sys.argv[1])
divs=[[] for _ in range(NMAX+1)]
for dd in range(1,NMAX+1):
    for k in range(dd,NMAX+1,dd): divs[k].append(dd)
res={}
for N in range(2,NMAX+1):
    D=divs[N]
    if len(D)<5: continue
    w=N**(-0.2)
    for i in range(len(D)):
        hi=D[i]*(1+w)
        j=i+1
        while j<len(D) and D[j]<=hi: j+=1
        cand=D[i+1:j]
        if len(cand)<4: continue
        for rest in itertools.combinations(cand,4):
            t=(D[i],)+rest
            g=0; L=1
            for v in t: g=gcd(g,v); L=L*v//gcd(L,v)
            if g!=1 or L!=N: continue
            G=math.log(N*t[4]/t[0])/(-math.log((t[4]-t[0])/t[0]))
            if G>5: continue
            k=int(math.floor(math.log2(N)))
            c,mg,arg=res.get(k,(0,99,None)); 
            res[k]=(c+1,min(mg,G),t if G<mg else arg)
for k in sorted(res): print("range 2^%d: count=%d min gamma=%.4f at %s"%(k,res[k][0],res[k][1],res[k][2]))
