#!/usr/bin/env python3
"""Referee's independent exact checker (written from scratch, not importing check_families.py).

Angles: every point z of a cluster is rotated by the Gaussian unit that brings it closest to a
reference ray; we then work with the actual lattice points u_j z_j (which lie on the same circle),
check they are pairwise distinct AS LATTICE POINTS, and compute the exact arc they span using
theta = 2*asin(|z-w|/(2R)) between the two extreme points (extremes found by exact cross products).
High precision via Decimal(100 digits) with argument-reduced series.
"""
from decimal import Decimal as D, getcontext
import itertools, math
getcontext().prec = 110

def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
def nrm(a): return a[0]*a[0]+a[1]*a[1]
def ggcd(a,b):
    while b != (0,0):
        n = nrm(b); t = mul(a,cj(b))
        q = ((2*t[0]+n)//(2*n), (2*t[1]+n)//(2*n))
        r = mul(b,q); a,b = b,(a[0]-r[0],a[1]-r[1])
    return a
def gdiv(a,b):
    n = nrm(b); t = mul(a,cj(b)); assert t[0]%n==0 and t[1]%n==0; return (t[0]//n,t[1]//n)
UNITS=[(1,0),(0,1),(-1,0),(0,-1)]

def dsqrt(x): return x.sqrt()
def dasin(x):
    # asin via atan: asin x = atan(x/sqrt(1-x^2)); atan by halving reduction
    return datan(x/dsqrt(1-x*x))
def datan(x):
    # reduce: atan x = 2 atan(x/(1+sqrt(1+x^2)))
    k=0
    while abs(x) > D("1e-6"):
        x = x/(1+dsqrt(1+x*x)); k+=1
    s=D(0); t=x; i=0; x2=x*x
    while True:
        term=t/(2*i+1)
        if term == 0 or abs(term) <= abs(s)*D(10)**-(getcontext().prec-5): break
        s += term if i%2==0 else -term
        t*=x2; i+=1
    return s*(2**k)

def cluster_C(points, reduce_gcd=True):
    """points: Gaussian integers of equal norm. Rotate each by a unit toward points[0],
    divide by the cluster gcd, verify distinct lattice points, return (C, N, k, gcdnorm, pts)."""
    N0 = nrm(points[0]); assert all(nrm(p)==N0 for p in points), "unequal norms"
    ref = points[0]
    rot=[]
    for p in points:
        # choose unit u maximizing Re(u p conj(ref))  (closest to ref direction)
        best=max(UNITS, key=lambda u: mul(mul(u,p),cj(ref))[0])
        rot.append(mul(best,p))
    g = rot[0]
    for p in rot[1:]: g = ggcd(g,p)
    if not reduce_gcd: g=(1,0)
    P=[gdiv(p,g) for p in rot]
    assert len(set(P))==len(P), "points coincide after unit rotation"
    N=nrm(P[0])
    # all within a half-plane around ref direction: Re(p conj(P0))>0
    for p in P: assert mul(p,cj(P[0]))[0] > 0
    # order by angle relative to P0 using cross product sign: key = Im/Re of p*conj(P0) (exact Fraction compare)
    from fractions import Fraction
    key=lambda p: Fraction(mul(p,cj(P[0]))[1], mul(p,cj(P[0]))[0])
    Ps=sorted(P,key=key)
    lo,hi=Ps[0],Ps[-1]
    R=dsqrt(D(N)); ch=dsqrt(D(nrm((hi[0]-lo[0],hi[1]-lo[1]))))
    theta=2*dasin(ch/(2*R))
    C=theta*dsqrt(R)
    return C,N,len(P),nrm(g),Ps

def best_k(points,k):
    C,N,kk,gn,Ps=cluster_C(points)
    R=dsqrt(D(N)); best=None
    for i in range(len(Ps)-k+1):
        lo,hi=Ps[i],Ps[i+k-1]
        ch=dsqrt(D(nrm((hi[0]-lo[0],hi[1]-lo[1]))))
        c=2*dasin(ch/(2*R))*dsqrt(R)
        best=c if best is None or c<best else best
    return best

# ---------------- golden family ----------------
fib=[0,1]
while len(fib)<2000: fib.append(fib[-1]+fib[-2])
def J(r): return (fib[r+1],fib[r])
F=(1,2)
def golden8(n):
    pts=[]
    for w in itertools.product([0,1],repeat=4):
        p=sum(w)
        if p%2==0: continue
        z = F if p==1 else cj(F)
        for s,ws in enumerate(w):
            b=J(n+s); z=mul(z, b if ws else cj(b))
        pts.append(z)
    return pts
s5=dsqrt(D(5)); phi=(1+s5)/2
C8=dsqrt(D(2))*(4+s5)/dsqrt(s5); C7=dsqrt(D(2))*phi**2*dsqrt(s5)
print("closed forms C8=%s C7=%s"%(str(C8)[:16],str(C7)[:16]))
for n in list(range(1,14))+[100,301]:
    C,N,k,gn,Ps=cluster_C(golden8(n))
    c7=best_k(golden8(n),7)
    print("golden n=%4d n%%3=%d k=%d gcdnorm=%d C8=%.12f (C8-lim=%.2e) C7=%.12f (C7-lim=%.2e) log10N=%.1f"%(
        n,n%3,k,gn,C,C-C8,c7,c7-C7,math.log10(N)))
# sanity: n=4 member has N=537606725
C,N,k,gn,Ps=cluster_C(golden8(4)); print("n=4: N=",N,"points",Ps)
