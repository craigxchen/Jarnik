"""Thm 3.6 numbers (single-polynomial barrier), T = 1 slab box, Sylvester-type lam = n/(n-1).
Conditions at Y^j (1<=|j|<=J): tau_j = floor(mu(1 + 4 ln2/W - |j| lam/sqrt n)) + 1, J = floor(1.07 sqrt(n)/lam);
count = sum_{1<=|j|<=J} binom(tau_j + d - 1, d)  (exact big integers);
available dimension >= |A| - binom(mu+d-1,d) >= ceil(3^d/2) - binom(mu+d-1,d)  (exact).
Reports, for mu = ceil(4w) and m_max, the margin, and mu_count = largest mu for which the barrier holds."""
import math
from mpmath import mp, mpf, sqrt, log, ceil, floor
mp.dps=40
def m_max_exact(d):
    lo,hi=0,2*d; t=3**d
    while lo<hi:
        mid=(lo+hi+1)//2
        if 8*math.comb(mid+d-1,d)<=t: lo=mid
        else: hi=mid-1
    return lo
def holds(n,mu,W):
    d=n-1; lam=mpf(n)/(n-1)
    J=int(floor(mpf('1.07')*sqrt(n)/lam))
    cnt=0
    for j in range(1,J+1):
        tau=int(floor(mu*(1+4*log(2)/W - j*lam/sqrt(n))))+1
        if tau<=0: break
        cnt+=2*math.comb(tau+d-1,d)
    avail=(3**d+1)//2 - math.comb(mu+d-1,d)
    return cnt<avail, cnt, avail
for n in [2966,3000,4096,8192]:
    d=n-1; W=3*n*log(n)-14*n
    w=sqrt(8*d*log(4*n)); m4=int(ceil(4*w)); mm=m_max_exact(d)
    out=[]
    for mu in (m4,mm):
        ok,cnt,av=holds(n,mu,W)
        out.append(f"mu={mu}: holds={ok}, ln(count)={math.log(cnt):.1f}, ln(avail)={math.log(av):.1f}")
    # mu_count: largest mu (search upward from m_max) with holds
    mu=mm
    while holds(n,mu+1,W)[0]: mu+=1
    print(f"n={n}: "+"; ".join(out)+f"; mu_count={mu} (rho(A) <= 2d = {2*d})")
