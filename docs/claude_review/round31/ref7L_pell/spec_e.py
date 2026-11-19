import sys; sys.path.insert(0,'.')
from pell_core import *
from run_b import build, A, beta
import math
from fractions import Fraction as Fr
def gi_divmod(a,b):  # Gaussian integers as (x,y) ints
    (x,y),(u,v)=a,b; n=u*u+v*v
    re=x*u+y*v; im=y*u-x*v
    qx=(2*re+n)//(2*n); qy=(2*im+n)//(2*n)
    return (qx,qy),(x-(qx*u-qy*v),y-(qx*v+qy*u))
def gi_gcd(a,b):
    while b!=(0,0): _,r=gi_divmod(a,b); a,b=b,r
    return a
def ev(p,Tv):
    acc=G(0)
    for c in reversed(p): acc=acc*G(Tv)+c
    return acc
def pext(p,q):  # returns (g,s,u) with s p + u q = g monic gcd
    r0,r1=trim(p),trim(q); s0,s1=P(1),[]; u0,u1=[],P(1)
    while r1:
        qq,rr=pdivmod(r0,r1)
        r0,r1=r1,rr; s0,s1=s1,psub(s0,pmul(qq,s1)); u0,u1=u1,psub(u0,pmul(qq,u1))
    inv=[r0[-1].inv()]
    return pmul(r0,inv),pmul(s0,inv),pmul(u0,inv)
def bezout(marks):
    g=marks[0]; coeffs=[P(1)]
    for z in marks[1:]:
        g2,s,u=pext(g,z); coeffs=[pmul(c,s) for c in coeffs]+[u]; g=g2
    return g,coeffs
def denom_lcm(polys):
    L=1
    for p in polys:
        for c in p: L=L*c.a.denominator//math.gcd(L,c.a.denominator); L=L*c.b.denominator//math.gcd(L,c.b.denominator)
    return L
for n,seeds in [(1,[Fr(1,2),Fr(1,3),Fr(2,7),Fr(-5,3)]),(1,[1,2,3,5]),(2,[Fr(1,2),Fr(1,3),Fr(2,7),Fr(-5,3)]),(2,[1,2,3,5]),(3,[1,2,3,5])]:
    m,_,_=build(A,beta,n,[Fr(x) for x in seeds],[(0,1),(0,2),(1,2),(0,3)])
    D=deg(m[0]); g,co=bezout(m)
    assert g==P(1)
    # check identity
    tot=[]
    for c,z in zip(co,m): tot=padd(tot,pmul(c,z))
    assert tot==P(1)
    Mz=denom_lcm(m); Mc=denom_lcm(co)
    # integer identity: sum (Mc*c)(Mz*z) = Mc*Mz ; so gcd of specialised Mz*z divides Mc*Mz
    print(f"n={n} seeds={seeds} D={D}  denominators: marks {Mz}, Bezout {Mc}; specialised gcd must divide {Mz*Mc}")
    for Tv in [10,100,1000,10**4,10**6,10**9]:
        vals=[]
        for z in m:
            w=ev(z,Tv); vals.append((int(w.a*Mz),int(w.b*Mz)))
            assert w.a*Mz==int(w.a*Mz) and w.b*Mz==int(w.b*Mz)
        gg=vals[0]
        for v in vals[1:]: gg=gi_gcd(gg,v)
        gn=gg[0]**2+gg[1]**2
        assert (Mz*Mc)**2 % gn ==0
        H=(vals[0][0]**2+vals[0][1]**2)//gn   # primitive norm R^2
        assert all((v[0]**2+v[1]**2)//gn==H for v in vals)
        R=math.sqrt(H) if H<2**1000 else None
        lR=0.5*math.log(H)
        arc=max(math.log(math.hypot(a[0]-b[0],a[1]-b[1]))-0.5*math.log(gn) for a,b in itertools.combinations(vals,2))
        print(f"   T={Tv:>10}: |gcd|^2={gn:>6}  log R={lR:9.3f}  exponent log(arc)/log R={arc/lR:.6f}  1/2+1/(2D)={0.5+0.5/D:.6f}  arc/sqrt(R T)={math.exp(arc-0.5*lR-0.5*math.log(Tv)):.4f}")
