# Independent exact computation of the minimal chord on |z|^2 = q_P, q_P = prod of primes p = r mod m in [P,2P].
import sys
from math import isqrt
def primes(lo,hi):
    return [p for p in range(max(lo,2),hi+1) if all(p%d for d in range(2,isqrt(p)+1))]
def gauss(p):  # a+bi with a^2+b^2=p
    for a in range(1,isqrt(p)+1):
        b2=p-a*a; b=isqrt(b2)
        if b*b==b2: return (a,b)
def mul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def toQ1(z):  # rotate by unit into x>0, y>=0
    x,y=z
    for _ in range(4):
        if x>0 and y>=0: return (x,y)
        x,y=-y,x
    raise
def run(P,mod,res):
    ps=[p for p in primes(P,2*P) if p%mod==res]
    q=1
    for p in ps: q*=p
    s=len(ps)
    pts=[(1,0)]
    for p in ps:
        a,b=gauss(p); n=[]
        for z in pts: n.append(mul(z,(a,b))); n.append(mul(z,(a,-b)))
        pts=n
    pts=[toQ1(z) for z in pts]
    assert all(x*x+y*y==q for x,y in pts)
    assert len(set(pts))==2**s
    K=4*q.bit_length()+20
    pts.sort(key=lambda z:(z[1]<<K)//z[0])
    # verify strict angular order exactly by cross products
    for i in range(len(pts)-1):
        (x1,y1),(x2,y2)=pts[i],pts[i+1]
        assert x1*y2-y1*x2>0
    best=None
    for i in range(len(pts)):
        z1=pts[i]
        z2=pts[i+1] if i+1<len(pts) else (-pts[0][1],pts[0][0])  # i*first
        d=(z1[0]-z2[0])**2+(z1[1]-z2[1])**2
        if best is None or d<best: best=d
    # floor(dmin/q^{1/4}) = largest k with k^4 q <= dmin^4 = best^2
    k=int((best**2/q)**0.25)+2
    while k**4*q>best**2: k-=1
    assert (k+1)**4*q>best**2
    print(P,mod,res,"s=",s,"orbits=",2**s,"dmin^2=",best,"floor=",k,"M=1 at C=8:",best**2>8**4*q)
for P in map(int,sys.argv[1:]):
    run(P,4,1)
