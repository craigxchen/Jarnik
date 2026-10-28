import sys; sys.path.insert(0,'.')
from gi import *
from math import isqrt, atan2, log, sqrt, comb, pi
from itertools import combinations
from fractions import Fraction as F
# (1) rho_{M,k} vs bound
bad=0
for k in range(2,9):
    for M in range(k+1,400):
        rho=min(F(comb(S,k)+comb(M-S,k),comb(M,k)) for S in range(1,M))
        if M>k:
            bd=F(2,2**k)*(1-F(k*(k-1),2*(M-k)))
            if rho<bd: bad+=1
        if M>=k+2*k*(k-1):
            # P(W0>=2^-k W) >= (rho-2^-k)/(1-2^-k) >= 2^-k-1 ?
            th=F(1,2**k)
            if (rho-th)/(1-th) < F(1,2**(k+1)): bad+=1
print('Lemma 3.2 numeric violations (k<=8, M<400):',bad)
print('rho_{10,5} =',float(min(F(comb(S,5)+comb(10-S,5),comb(10,5)) for S in range(1,10))))
# (2) the 10-point cluster on n=1176852625
n=1176852625; R=sqrt(sqrt(n))**2  # R = sqrt(n)
R=sqrt(n)
P=[]
for x in range(0,isqrt(n)+1):
    y2=n-x*x; y=isqrt(y2)
    if y*y==y2:
        for sx in (1,-1):
            for sy in (1,-1): P.append((sx*x,sy*y))
P=sorted(set(P),key=lambda z: atan2(z[1],z[0]))
L=len(P); best=None
for i in range(L):
    W=[P[(i+j)%L] for j in range(10)]
    d=max(sqrt(gnorm((a[0]-b[0],a[1]-b[1]))) for a,b in combinations(W,2))
    if best is None or d<best[0]: best=(d,W)
d,W=best; C=d/sqrt(R)
print('points on circle',L,' best 10-window diameter/sqrt(R) =',C)
g=ggcdlist(W); print('cluster gcd norm',gnorm(g))
Wt=log(n); vals=[]; minexp=9; maxexp=0; ok=True
for T in combinations(W,5):
    gg=ggcdlist(list(T)); W0=log(gnorm(gg)); vals.append(W0/Wt)
    Rpp=R/sqrt(gnorm(gg))
    if Rpp < sqrt(R)-1e-9: ok=False
    zeta=[gexact(z,gg) for z in T]
    dd=max(sqrt(gnorm((a[0]-b[0],a[1]-b[1]))) for a,b in combinations(zeta,2))
    e=log(dd)/log(Rpp); minexp=min(minexp,e); maxexp=max(maxexp,e)
print('E[W0/W] over 5-subsets =',sum(vals)/len(vals),' Rpp>=sqrt(R) always:',ok,' reduced exponents range',round(minexp,3),round(maxexp,3))
