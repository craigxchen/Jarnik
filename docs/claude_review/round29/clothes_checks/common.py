import math, sympy
from math import gcd, isqrt
def gpoints(N):
    """all z=(x,y) with x^2+y^2=N, via Gaussian factorisation (N odd, primes 1 mod 4 or squares of 3 mod 4)."""
    f = sympy.factorint(N)
    pts=[(1,0)]
    def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    for p,e in f.items():
        if p==2:
            new=[]
            for z in pts:
                w=z
                for _ in range(e): w=mul(w,(1,1))
                new.append(w)
            pts=new; continue
        if p%4==3:
            if e%2: return []
            pts=[(z[0]*p**(e//2), z[1]*p**(e//2)) for z in pts]; continue
        # find pi
        a=next(a for a in range(1,isqrt(p)+1) if isqrt(p-a*a)**2==p-a*a)
        b=isqrt(p-a*a); pi=(a,b); pib=(a,-b)
        pw=[(1,0)]; pwb=[(1,0)]
        for _ in range(e): pw.append(mul(pw[-1],pi)); pwb.append(mul(pwb[-1],pib))
        pts=[mul(z,mul(pw[k],pwb[e-k])) for z in pts for k in range(e+1)]
    units=[(1,0),(0,1),(-1,0),(0,-1)]
    S=set(mul(z,u) for z in pts for u in units)
    assert all(x*x+y*y==N for x,y in S)
    return sorted(S)
def arcmax(N,C):
    """max # points on an arc of length <= C*N^{1/4} (chord criterion: angular width C N^{-1/4})"""
    P=gpoints(N); R=math.sqrt(N); w=C*N**0.25/R
    ang=sorted(math.atan2(y,x)%(2*math.pi) for x,y in P)
    n=len(ang); ang2=ang+[a+2*math.pi for a in ang]; best=0; j=0
    for i in range(n):
        while j<i+n and ang2[j]-ang2[i]<=w+1e-15: j+=1
        best=max(best,j-i)
    return best
