# Independent census of D_4(N) over odd split-supported N <= X.
import sys, math, numpy as np
from math import isqrt, gcd, log
from itertools import combinations
X=int(sys.argv[1]); B=4
spf=np.zeros(X+1,dtype=np.int32)
for p in range(2,isqrt(X)+1):
    if spf[p]==0:
        sl=spf[p*p::p]; m=(sl==0); sl[m]=p; spf[p*p::p]=sl
def split_supported(n):
    if n%2==0: return False
    while n>1:
        p=int(spf[n]) or n
        if p%4!=1: return False
        while n%p==0: n//=p
    return True
cand={}
nviol_upper=0
for s in range(1,isqrt(X)+1):
    for t in range(1,s):
        d=s*s+t*t
        if d>X: break
        if 256*t**4>d: break
        if (s+t)%2==0 or gcd(s,t)!=1: continue
        assert isqrt(d)==s
        assert split_supported(d)
        mmax=min(X//d, d//(256*t**4))
        for m in range(1,mmax+1):
            if m%2==0 or not split_supported(m): continue
            N=d*m
            # exact conditions
            assert 256*t**4*N<=d*d and 256*N<=d*d
            if not (d*d <= N*math.log(N)**(2*B)):  # d <= sqrt N W^B
                nviol_upper+=1; continue
            cand.setdefault(N,[]).append((d,s,t))
from collections import Counter
cnt=Counter(len(v) for v in cand.values())
print("X",X,"parents by #candidates",sorted(cnt.items()),"upper-cutoff rejections",nviol_upper)
mx=max(cnt); print("max parents",[ (N,v) for N,v in cand.items() if len(v)==mx][:5])
tri=0
for N,v in cand.items():
    for a,b,c in combinations(v,3):
        P=a[0]*b[0]*c[0]
        if isqrt(P)**2==P: tri+=1
print("square-product triples",tri)
