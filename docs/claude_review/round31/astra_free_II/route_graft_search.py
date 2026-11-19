"""Optional floating-point discovery search. Use graft_check.py for certification."""
import math
import sys
from fractions import Fraction
import numpy as np

def candidates(B, p, j, phi):
    rows=[]
    for a in range(1, math.isqrt(B)+1):
        b=np.arange(0,math.isqrt(B-a*a)+1,dtype=np.int64)
        b=b[(np.gcd(a,b)==1)&((a*a+b*b)%p!=0)]
        if len(b):
            rows.append(np.column_stack((np.full(len(b),a),b)))
    z=np.concatenate(rows)
    angles=(np.arctan2(z[:,1],z[:,0])+j*phi)%(math.pi/2)
    perm=np.argsort(angles)
    return angles[perm], z[perm]

def search(X, p):
    phi=math.atan(math.sqrt(2))
    t1,z1=candidates(X//p,p,1,phi)
    t2,z2=candidates(X//(p*p),p,2,phi)
    k=np.searchsorted(t1,t2)
    best=None
    for offset in (-1,0):
        ix=np.clip(k+offset,0,len(t1)-1)
        delta=t1[ix]-t2
        for a in np.where(abs(delta)<0.2/X)[0]:
            h1,h2=t1[ix[a]],t2[a]
            mid=(h1+h2)/2
            swapped=mid>math.pi/4
            t=math.tan(math.pi/2-mid if swapped else mid)
            rr=Fraction(t).limit_denominator(math.floor(math.sqrt(X/(1+t*t))))
            v0=(rr.numerator,rr.denominator) if swapped else (rr.denominator,rr.numerator)
            x0=sum(u*u for u in v0)
            if x0>X or x0%p==0: continue
            h0=math.atan2(v0[1],v0[0])
            v1,v2=z1[ix[a]],z2[a]
            heights=(x0,p*(v1@v1),p*p*(v2@v2))
            diameter=max(heights)*(max(h0,h1,h2)-min(h0,h1,h2))
            if best is None or diameter<best[0]:
                best=(diameter,tuple(v0),tuple(v1),tuple(v2),tuple(heights),(h0,h1,h2))
    print(best)

if __name__=='__main__':
    search(int(sys.argv[1]),int(sys.argv[2]))
