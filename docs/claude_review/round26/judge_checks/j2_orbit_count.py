"""Judge check J2 (exact big integers; J and tau via exact-ish mpmath at 50 digits, same conventions as ze.md Thm 3.6).
Recompute the referee's AGL-orbit-inclusive count at n = 4096 (k=12): cutoff mu where
|AGL(12,2)| * sum_{1<=|j|<=J} binom(tau_j+d-1,d) < ceil(3^d/2) - binom(mu+d-1,d) stops holding; and the powers-only count."""
import math
from mpmath import mp, mpf, sqrt, log, floor
mp.dps=50
n=4096; d=n-1; k=12
W=3*n*log(n)-14*n; lam=mpf(n)/(n-1)
J=int(floor(mpf('1.07')*sqrt(n)/lam))
agl=2**k
for i in range(k): agl*=(2**k-2**i)
T3=3**d
def counts(mu):
    cnt=0
    for j in range(1,J+1):
        tau=int(floor(mu*(1+4*log(2)/W-j*lam/sqrt(n))))+1
        if tau<=0: break
        cnt+=2*math.comb(tau+d-1,d)
    avail=(T3+1)//2-math.comb(mu+d-1,d)
    return cnt,avail
def mmax():
    lo,hi=0,2*d
    while lo<hi:
        mid=(lo+hi+1)//2
        if 8*math.comb(mid+d-1,d)<=T3: lo=mid
        else: hi=mid-1
    return lo
mm=mmax(); print("J =",J,"m_max =",mm)
for mu in (2256,2558,2559,2560,mm):
    c,a=counts(mu)
    print(f"mu={mu}: powers-only holds={c<a}; orbit holds={agl*c<a}; ln(agl*c)={math.log(agl*c):.2f} ln(avail)={math.log(a):.2f}")
