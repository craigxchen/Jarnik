"""Evidence for Heuristic 6.3: for Haar-random phi in T^s (T = R/(pi/2)Z), estimate the probability that ALL monomial
gap events G_v (overview Lemma 5.1) hold for every v (up to sign) with w(v) <= WMAX, and compare with the Poisson
prediction exp(-E), E = sum_{v up to sign, w(v)<=WMAX} Prob(Bad_v). Also reports the tail intensity beyond WMAX
(upper bound via the product formula) and, as a control, whether the genuine angles theta_p pass (they must).
Usage: python3 poisson_gaps.py WMAX NSAMP p1 p2 ...
"""
import sys, math
import numpy as np
from glib import gprime
WMAX=float(sys.argv[1]); NS=int(sys.argv[2]); ps=[int(x) for x in sys.argv[3:]]
s=len(ps); lw=[math.log(p) for p in ps]
vs=[]
def dfs(j,v,w,started):
    if j==s:
        if started: vs.append((tuple(v),w))
        return
    m=0
    while w+abs(m)*lw[j]<=WMAX+1e-12:
        for sg in ([1] if m==0 else ([1] if not started else [1,-1])):
            v.append(sg*m); dfs(j+1,v,w+abs(m)*lw[j],started or m!=0); v.pop()
        m+=1
dfs(0,[],0.0,False)
V=np.array([v for v,w in vs],dtype=np.float64); Wv=np.array([w for v,w in vs])
a=np.arcsin(np.exp(-Wv/2)); b=np.arcsin(np.exp(-Wv/2)/math.sqrt(2))
pbad=(4/math.pi)*(a+b)
E=pbad.sum()
# tail bound: (1/2) sum_{v!=0} 2.18 e^{-w/2} minus enumerated part
full=0.5*(np.prod([1+2/(math.sqrt(p)-1) for p in ps])-1)
tail=2.18*(full-0.5*0) - 2.18*np.exp(-Wv/2).sum()
Q=math.pi/2
def check(phi, tol=0.0):  # phi: (n,s); tol>0 only for the genuine control (equality cases |Im A|=1)
    t=np.mod(phi@V.T, Q)                    # (n, nv)
    d0=np.minimum(t, Q-t); d1=np.abs(t-Q/2)  # dist to (pi/2)Z and to pi/4+(pi/2)Z
    bad=(d0<a[None,:]*(1-tol))|(d1<b[None,:]*(1-tol))
    return ~bad.any(axis=1)
rng=np.random.default_rng(12345)
good=0; done=0; CH=max(1,int(2e7//max(1,len(vs))))
while done<NS:
    n=min(CH,NS-done); phi=rng.random((n,s))*Q
    good+=int(check(phi).sum()); done+=n
th=np.array([[math.atan2(gprime(p)[1],gprime(p)[0]) for p in ps]])
print("primes",ps,"s=%d sum p^-1/2=%.3f"%(s,sum(p**-0.5 for p in ps)))
print("WMAX=%.1f #v(up to sign)=%d E_enum=%.4f tail_intensity<=%.4f"%(WMAX,len(vs),E,max(tail,0)))
print("Haar: P(all enumerated gaps hold)=%.5f (+-%.5f)  Poisson exp(-E)=%.5f  ratio=%.3f"%(good/NS, math.sqrt(good)/NS, math.exp(-E), (good/NS)/math.exp(-E)))
print("genuine angles pass (rel. tol 1e-9 for equality cases):", bool(check(th,1e-9)[0]))
