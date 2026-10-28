# Referee: independent recomputation of prove.md §7.3 (#v, E, direct Monte Carlo of Haar measure of Good(WMAX)),
# plus Prop 6.1(2) product formula and the genuine-angle gap check.
import math, itertools, sys
import numpy as np
from mpmath import mp, mpf, atan2, asin, exp, sqrt, pi as mpi, nint
mp.dps=40
pi=math.pi
def split_primes_from(a, cnt):
    out=[]; p=a
    while len(out)<cnt:
        if p%4==1 and all(p%d for d in range(2,int(p**.5)+1)): out.append(p)
        p+=1
    return out
def enum_v(primes, WMAX):
    L=[math.log(p) for p in primes]; s=len(primes); out=[]
    def rec(j, cur, w):
        if j==s:
            if any(cur): out.append((tuple(cur), w))
            return
        m=0
        while w+abs(m)*L[j]<=WMAX+1e-12:
            for sg in ([m] if m==0 else [m,-m]):
                cur.append(sg); rec(j+1,cur,w+abs(m)*L[j]); cur.pop()
            m+=1
    rec(0,[],0.0)
    # up to sign: keep v whose first nonzero coordinate is positive
    out=[(v,w) for (v,w) in out if next(x for x in v if x!=0)>0]
    return out
def pbad(w):
    x=math.exp(-w/2); return (4/pi)*(math.asin(x)+math.asin(x/math.sqrt(2)))
def gauss(p):
    for a in range(1,int(p**.5)+1):
        b2=p-a*a; b=int(round(b2**.5))
        if b*b==b2: return a,b
def run(primes, WMAX, nsamp, seed=0):
    V=enum_v(primes,WMAX); V.sort(key=lambda t:t[1])
    E=sum(pbad(w) for _,w in V)
    sp=sum(p**-0.5 for p in primes)
    prodf=1.0
    for p in primes: prodf*=1+2/(math.sqrt(p)-1)
    # genuine angles check with mpmath
    th=[atan2(gauss(p)[1],gauss(p)[0]) for p in primes]
    worst=mpf(10)
    for v,w in V:
        t=sum(vi*ti for vi,ti in zip(v,th))
        d0=abs(t-(mpi/2)*nint(t/(mpi/2))); u=t-mpi/4; d1=abs(u-(mpi/2)*nint(u/(mpi/2)))
        r=min(d0/asin(exp(-mpf(w)/2)), d1/asin(exp(-mpf(w)/2)/sqrt(2)))
        worst=min(worst,r)
    # direct MC with early elimination
    rng=np.random.default_rng(seed); alive_total=0; chunk=200000
    Vm=np.array([v for v,_ in V],dtype=np.float64); a0=np.array([math.asin(math.exp(-w/2)) for _,w in V]); a1=np.array([math.asin(math.exp(-w/2)/math.sqrt(2)) for _,w in V])
    done=0
    while done<nsamp:
        n=min(chunk,nsamp-done); phi=rng.random((n,len(primes)))*(pi/2); alive=np.ones(n,bool)
        for i in range(len(V)):
            idx=np.nonzero(alive)[0]
            if idx.size==0: break
            t=(phi[idx]@Vm[i])%(pi/2)
            d0=np.minimum(t,pi/2-t); u=(t-pi/4)%(pi/2); d1=np.minimum(u,pi/2-u)
            bad=(d0<a0[i])|(d1<a1[i]); alive[idx[bad]]=False
        alive_total+=int(alive.sum()); done+=n
    P=alive_total/nsamp
    se=math.sqrt(max(P*(1-P),1e-300)/nsamp)
    print(f"primes={primes} sum p^-1/2={sp:.3f} WMAX={WMAX} #v={len(V)} E={E:.3f} 1.08*(prod-1)={1.08*(prodf-1):.3f} "
          f"(1/2)*2.1736*(prod-1)={1.0868*(prodf-1):.3f} | MC n={nsamp} good={alive_total} P={P:.3e} logP={math.log(P) if P>0 else float('nan'):.3f} "
          f"(+-{se/P if P>0 else 0:.2%}) -logP/E={(-math.log(P)/E) if P>0 else float('nan'):.3f} | genuine worst gap ratio={float(worst):.6f}", flush=True)
run([101,109,113,137,149,157],24,400000,1)
run([53,61,73],22,400000,2)
run(split_primes_from(53,5),22,1000000,3)
run(split_primes_from(53,7),22,2000000,4)
run([29,37,41,53,61,73,89,97],22,8000000,5)
