"""Item (4c): what (IX.1) gives with the best Gaussian direction b for a generic arc direction.
F(theta) = min_b |b| (1 + sqrt(R') |sin(theta - arg b)|) over nonzero b in Z[i].
X_b <<_C F(theta), and (IX.1) needs F = o(s^{3/2}), s <= (1+o(1)) log N'/loglog N'.
We compute F by continued fractions (convergents and semiconvergents of the slope), checked by
brute force at small R'."""
import math, random
import numpy as np
from mpmath import mp, mpf, atan2, sin, sqrt, cos, tan, pi
mp.dps=60
def F_brute(theta,Rp,B):
    a=np.arange(-B,B+1); A,Bb=np.meshgrid(a,a); A=A.ravel(); Bb=Bb.ravel()
    m=(A*A+Bb*Bb>0)&(A*A+Bb*Bb<=B*B)
    A=A[m].astype(float); Bb=Bb[m].astype(float)
    nb=np.hypot(A,Bb); s=np.abs(np.sin(float(theta)-np.arctan2(Bb,A)))
    return float(np.min(nb*(1+math.sqrt(Rp)*s)))
def F_cf(theta,Rp):
    # reduce to slope t = tan(theta mod pi/2) in [0,1] after symmetries; directions (q,p) with p/q ~ t
    th=theta % (pi/2)
    if th>pi/4: th=pi/2-th
    t=tan(th); best=mpf(10)**40
    # continued fraction of t with semiconvergents
    p0,q0,p1,q1=0,1,1,0; x=t
    for _ in range(80):
        a=int(mp.floor(x))
        for k in range(1,a+1):
            p,q=k*p1+p0,k*q1+q0
            nb=sqrt(p*p+q*q); val=nb*(1+sqrt(Rp)*abs(sin(th-atan2(p,q))))
            best=min(best,val)
        p0,q0,p1,q1=p1,q1,a*p1+p0,a*q1+q0
        nb=sqrt(p1*p1+q1*q1); best=min(best,nb*(1+sqrt(Rp)*abs(sin(th-atan2(p1,q1)))))
        if x==a: break
        x=1/(x-a)
        if q1>10*mp.sqrt(sqrt(Rp))*10: break
    return float(best)
random.seed(1)
# brute-force cross-check
for Rp in [1e6,1e8]:
    for _ in range(5):
        th=mpf(random.random())*2*pi
        fb=F_brute(th,Rp,int(4*Rp**0.25)+5); fc=F_cf(th,Rp)
        assert abs(fb-fc)<1e-6*fb, (fb,fc)
print("continued-fraction minimiser agrees with brute force on 10 samples")
for e in [8,12,16,24,40]:
    Rp=mpf(10)**e
    vals=sorted(F_cf(mpf(random.random())*2*pi,Rp) for _ in range(400))
    W=float(2*mp.log(Rp)); s=W/math.log(W)
    print(f"R'=1e{e}: R'^(1/4)={float(Rp**0.25):.3g}; F/R'^(1/4) quantiles 1%,50%,99%:",
          [round(vals[int(f*len(vals))]/float(Rp**0.25),3) for f in (0.01,0.5,0.99)],
          f"; s^(3/2) bound ~ {s**1.5:.3g}")
# proportion of directions with F <= T: heuristic count sum_{|b|<=T} 2T/(|b| sqrt R') / (2 pi), checked by sampling
Rp=mpf(10)**12; T=51.0
hits=sum(F_cf(mpf(random.random())*2*pi,Rp)<=T for _ in range(100000))
pred=0.0
for a in range(-int(T),int(T)+1):
    for b in range(-int(T),int(T)+1):
        nb=math.hypot(a,b)
        if 0<nb<=T and math.gcd(a,b)==1:
            pred+=2*(T/nb-1)/math.sqrt(float(Rp)) if T/nb>1 else 0
print(f"R'=1e12, T={T}: sampled proportion of theta with F<=T: {hits/100000:.4g}; first-order prediction {pred/(2*math.pi):.4g}")
