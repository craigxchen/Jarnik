"""Rebuild the primitive 5-tuple of a cluster (L; X_1..X_4) and print its cut profile."""
import math, itertools, sys
from glib import *
from sympy import factorint
from enum_chart import beta_bar
def tuple_of(L,Xs):
    Z=(1,0)
    for X in Xs: Z=lcm(Z,beta_bar(X,L))
    pts=[Z]
    for X in Xs:
        w=(X,L); pts.append(exdiv(mul(Z,w),conj(w)))
    H=nrm(Z); assert all(nrm(p)==H for p in pts)
    g=pts[0]
    for p in pts[1:]: g=gcd(g,p)
    assert nrm(g)==1
    d=max(abs(complex(*a)-complex(*b)) for a,b in itertools.combinations(pts,2))
    assert abs(d/math.sqrt(H)-2*L/math.sqrt(Xs[0]**2+L*L))<1e-9
    return pts,H
def profile(pts,H):
    W=math.log(H); fac=factorint(H); cuts={}
    for p,e in fac.items():
        pi=split_prime(p); a=[val(x,pi) for x in pts]
        for t in range(1,e+1):
            hi=frozenset(k for k in range(5) if a[k]>=t)
            side=hi if len(hi)<=2 else frozenset(range(5))-hi
            cuts[side]=cuts.get(side,0)+math.log(p)
    pw={e:cuts.get(frozenset(e),0)/W for e in itertools.combinations(range(5),2)}
    ws=sum(w for c,w in cuts.items() if len(c)==1)/W
    return pw,ws,fac
if __name__=='__main__':
    L=int(sys.argv[1]); Xs=[int(x) for x in sys.argv[2].split(',')]
    pts,H=tuple_of(L,Xs); pw,ws,fac=profile(pts,H)
    print('H=',H,'=',fac); print('points',pts)
    print('pair weights/W',{k:round(v,3) for k,v in pw.items()}); print('singleton/W',round(ws,3))
