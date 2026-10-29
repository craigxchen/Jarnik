# Pell points of mnn31, gamma -> 7/2, cone membership, m-czech not on cones, inversion images.
import math, itertools, sympy as sp
from math import gcd
from functools import reduce
def lcm(x,y): return x//gcd(x,y)*y
# (2u+5)^2 - 3(2w+3)^2 = -2 : X^2-3Y^2=-2, X=2u+5, Y=2w+3. fundamental: X=1,Y=1; unit 2+sqrt3
sols=[]; X,Y=1,1
for i in range(60):
    if X%2==1 and Y%2==1:
        uu=(X-5)//2; ww=(Y-3)//2
        if uu>0 and ww>0: sols.append((uu,ww))
    X,Y=2*X+3*Y, X+2*Y
def cone_val(v,base,pairs):
    (p,q),(r,s)=pairs
    return (v[p]-v[base])*(v[q]-v[base])-(v[r]-v[base])*(v[s]-v[base])
def cones():
    out=[]
    for base in range(5):
        rest=[i for i in range(5) if i!=base]
        r0=rest[0]
        for j in rest[1:]:
            p1=(r0,j); p2=tuple(x for x in rest if x not in p1)
            out.append((base,(p1,p2)))
    return out
CONES=cones(); assert len(CONES)==15
print("Pell points:",len(sols))
for (uu,ww) in sols[:25]:
    n=uu*(uu+2)*(uu+3)*(uu+5)
    Dv=[(uu+2)*(uu+3),(uu+2)*(uu+5),(uu+3)*(uu+5),3*(ww+1)*(ww+3),3*(ww+2)*(ww+3)]
    assert all(n%x==0 for x in Dv) and len(set(Dv))==5
    assert (uu+2)*(uu+3)==3*(ww+1)*(ww+2)
    on=[i for i,(bb,pp) in enumerate(CONES) if cone_val(Dv,bb,pp)==0]
    g=reduce(gcd,Dv); P=[x//g for x in Dv]; L=reduce(lcm,P)
    H=max(P)/min(P)*L; dist=(max(P)-min(P))/min(P)
    gam=math.log(H)/(-math.log(dist))
    # co-divisors (inversion): L/P_a
    Q=[L//x for x in P]
    on_inv=[i for i,(bb,pp) in enumerate(CONES) if cone_val(Q,bb,pp)==0]
    gQ=reduce(gcd,Q); LQ=reduce(lcm,Q)
    sq=math.isqrt(n)
    Cw=(max(Dv)-sq)/n**0.25
    print(f"u~1e{math.log10(uu):.1f} cones={on} gamma={gam:.5f} log(dist^3.5 H)={3.5*math.log(dist)+math.log(H):.3f} logL/logn={math.log(L)/math.log(n):.4f} C={Cw:.4f} inv:gcd={gQ} lcm_eq={LQ==L} inv_on_cones={on_inv}")
# m-czech (round23 Thm 2) generic W-member not on any cone: symbolic in x on conic, (a,b,c)=(1,2,4)
x,z=sp.symbols('x z')
def wcurve(A,Bb,Cc):
    Dd=A+Bb+Cc; r2=[A**2,Bb**2,Cc**2,Dd**2]
    S=sum(r2); e2=sum(p*q for p,q in itertools.combinations(r2,2)); e3=sum(p*q*s for p,q,s in itertools.combinations(r2,3))
    m=A*Bb*Cc*Dd; Bc=sp.Rational(e3-S*m,4)
    F=[(x-A)*(x-Bb)*(x-Cc)*(x+Dd), x**4-sp.Rational(S,2)*x**2+m+2*x*z,(x**2-Bb**2)*(x**2-Cc**2),(x**2-A**2)*(x**2-Cc**2),(x**2-A**2)*(x**2-Bb**2)]
    return F,m,Bc
F,m,Bc=wcurve(1,2,4); print("m-czech m,B:",m,Bc)
# check on Pell points z^2 = 56x^2+65 numerically-exactly
pts=[]
for xx in range(1,200000):
    zz2=m*xx*xx+Bc
    zz=math.isqrt(int(zz2))
    if zz*zz==zz2: pts.append((xx,zz))
print("m-czech integer points x<2e5:",pts[:8])
for (xx,zz) in pts[:6]:
    for sgn in (1,-1):
        Fv=[int(f.subs({x:xx,z:sgn*zz})) for f in F]
        on=[i for i,(bb,pp) in enumerate(CONES) if cone_val(Fv,bb,pp)==0]
        print("  x=",xx,"sign",sgn,"on cones:",on)
# symbolic: cone polys restricted to the W-curve (generic a,b,c random rationals), reduce mod conic
for (A,Bb,Cc) in [(1,2,4),(sp.Rational(3,2),sp.Rational(5,3),sp.Rational(7,5))]:
    F,m,Bc=wcurve(A,Bb,Cc)
    bad=0
    for (bb,pp) in CONES:
        e=sp.expand(cone_val(F,bb,pp))
        e=sp.expand(e.subs(z**2,m*x**2+Bc)); e=sp.expand(e.subs(z**2,m*x**2+Bc))
        p0=sp.Poly(e,x,z)
        if p0.is_zero: bad+=1
    print("W-curve",(A,Bb,Cc),"contained in",bad,"of 15 cones")
# inversion image of the cone: (u1-u2)(u1-u3)u4u5-(u1-u4)(u1-u5)u2u3, factor
u1,u2,u3,u4,u5=sp.symbols('u1:6')
inv=sp.factor((u1-u2)*(u1-u3)*u4*u5-(u1-u4)*(u1-u5)*u2*u3)
print("iota(cone) =",inv)
