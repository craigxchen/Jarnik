# Independent check of gauss.md Thm 4.1: instance X=-3,m=3 and Jacobian rank of (X,m,s)->T_4
import sympy as sp
from fractions import Fraction as F
X,m,s=sp.symbols('X m s')
def fourth(roots):
    s1=sum(roots); s2=roots[0]*roots[1]+roots[0]*roots[2]+roots[1]*roots[2]
    s1=sp.expand(s1); s2=sp.expand(s2)
    re1,im1=sp.re(s1),sp.im(s1); im2=sp.im(s2)
    return sp.simplify((re1*im1-im2)/im1) - sp.I*im1
def family(Xv,mv):
    xa=Xv; xb=Xv*(1-Xv)/(mv**2+Xv); xc=1-xa-xb; u=mv*xb
    r123=sp.I; r12=u+sp.I*(xa-1); r13=u+sp.I*(xb-1); r23=u+sp.I*(xc-1)
    r1=fourth([r12,r13,r123]); r2=fourth([r12,r23,r123]); r3=fourth([r13,r23,r123])
    return dict(r123=r123,r12=r12,r13=r13,r23=r23,r1=r1,r2=r2,r3=r3)
R=family(sp.Integer(-3),sp.Integer(3))
print({k:sp.nsimplify(v) for k,v in R.items()})
blocks={1:['r1','r12','r13','r123'],2:['r2','r12','r23','r123'],3:['r3','r13','r23','r123']}
for x,bl in blocks.items():
    P=sp.expand(sp.prod([s-R[b] for b in bl]))
    print('x=%d  Im P_x(s) ='%x, sp.simplify(sp.im(P.subs(s,sp.Symbol('t',real=True)))))
# Jacobian rank of t_x = Im P_x / Re P_x in (X,m,s) at rational points
Rs=family(X,m)
t=sp.Symbol('t',real=True)
coords=[]
Xr,mr=sp.symbols('Xr mr',real=True)
Rr=family(Xr,mr)
for x,bl in blocks.items():
    P=sp.expand(sp.prod([t-Rr[b] for b in bl]))
    coords.append(sp.re(P)/sp.im(P))   # Re/Im (Im P is s-independent constant on the family)
J=sp.Matrix([[sp.diff(c,v) for v in (Xr,mr,t)] for c in coords])
for pt in [(-3,3,10),(-3,3,1000),(2,5,7),(F(1,2),F(3,2),11),(-5,1,3),(7,2,-4)]:
    sub={Xr:sp.Rational(pt[0]),mr:sp.Rational(pt[1]),t:sp.Rational(pt[2])}
    Jn=J.subs(sub)
    print(pt,'rank',Jn.rank(), 'det', sp.nsimplify(Jn.det()))
