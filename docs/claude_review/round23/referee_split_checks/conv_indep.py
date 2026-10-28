"""Rational points of C_{abc} at convergents p/q of sqrt(m) (parameter t = slope through (a,z_a)):
exact primitive tuple, H=(max/min)*lcm, gamma, log(dist^{15/4}H). Independent of split_checks/w_points.py."""
from fractions import Fraction as Fr
from math import gcd, log, isqrt
import sys
def lcm(p,q): return p//gcd(p,q)*q
def run(A,Bv,Cv,ncv=70):
    a,b,c=Fr(A),Fr(Bv),Fr(Cv); d=a+b+c; sq=[a*a,b*b,c*c,d*d]; S=sum(sq)
    e3=sq[0]*sq[1]*sq[2]+sq[0]*sq[1]*sq[3]+sq[0]*sq[2]*sq[3]+sq[1]*sq[2]*sq[3]
    m=a*b*c*d; B=(e3-S*m)/4
    za=(a*a*b+a*a*c+a*b*b+a*c*c-b*b*c-b*c*c)/2
    assert za*za==m*a*a+B
    M=int(m); r=isqrt(M)
    # continued fraction of sqrt(M)
    p0,q0,p1,q1=1,0,r,1; mm,dd,aa=0,1,r; out=[]
    for k in range(ncv):
        t=Fr(p1,q1)
        x=(a*(t*t+m)-2*za*t)/(t*t-m); z=za+t*(x-a)
        assert z*z==m*x*x+B
        F=[(x-a)*(x-b)*(x-c)*(x+d), x**4-S/2*x*x+m+2*x*z,(x*x-b*b)*(x*x-c*c),(x*x-a*a)*(x*x-c*c),(x*x-a*a)*(x*x-b*b)]
        if any(f==0 for f in F): pass
        else:
            den=1
            for f in F: den=lcm(den,f.denominator)
            u=[int(f*den) for f in F]
            if all(v>0 for v in u) or all(v<0 for v in u):
                u=sorted(abs(v) for v in u); g=0
                for v in u: g=gcd(g,v)
                u=[v//g for v in u]; L=1
                for v in u: L=lcm(L,v)
                if len(set(u))==5:
                    lH=log(L)+log(u[-1])-log(u[0]); ld=log(u[-1]-u[0])-log(u[0])
                    out.append((lH,lH/(-ld),3.75*ld+lH))
        mm=dd*aa-mm; dd=(M-mm*mm)//dd; aa=(r+mm)//dd
        p0,q0,p1,q1=p1,q1,aa*p1+p0,aa*q1+q0
    return out
for abc in [(1,2,4),(1,3,5),(2,3,7),(1,5,6)]:
    o=run(*abc)
    f=[r for r in o if r[1]<4]
    print(abc,"pts",len(o),"first gamma<4 at logH=%.1f"%f[0][0] if f else "none",
          "logH max=%.0f"%max(r[0] for r in o),"range log(dist^3.75 H)=[%.2f,%.2f]"%(min(r[2] for r in o[5:]),max(r[2] for r in o[5:])),
          "last gamma=%.4f"%o[-1][1])
