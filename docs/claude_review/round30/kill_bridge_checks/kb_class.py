# Generic independent C_inf(y,eps) for the full 8-point cube cluster (upper bound: sampled local maxima).
import sys, itertools, math, sympy as sp
from fractions import Fraction as Fr
import kb_gauss as KG
from kb_local import badprimes, vl
def gprimes(l):
    if l==2: return [(1,1)]
    if l%4==3: return [(l,0)]
    for x in range(1,math.isqrt(l)+1):
        s=math.isqrt(l-x*x)
        if s*s==l-x*x: return [(x,s),(x,-s)]
def vpi_gint(z,pi):
    n=pi[0]**2+pi[1]**2; v=0
    while True:
        x=z[0]*pi[0]+z[1]*pi[1]; yy=z[1]*pi[0]-z[0]*pi[1]
        if x%n or yy%n: return v
        z=(x//n,yy//n); v+=1
def vpi(val,pi,l):
    a,b=val; den=a.denominator*b.denominator//math.gcd(a.denominator,b.denominator)
    return vpi_gint((int(a*den),int(b*den)),pi)-vl(Fr(den),l)*(2 if l==2 else 1)
def logg(rts,t0,l):
    vals=KG.values(rts,t0)
    return sum(math.log(p[0]**2+p[1]**2)*min(vpi(v,p,l) for v in vals) for p in gprimes(l))
def localmax(rts,l):
    dens=set()
    for a,b in rts: dens|={a.denominator,b.denominator}
    jl=max(vl(Fr(d),l) for d in dens)
    if l in (2,3,5) or jl>0:
        jm=jl+2; k={2:8,3:5,5:4}.get(l,3)
        S=[Fr(a,l**j) for j in range(jm+1) for a in range(l**(j+k)) if j==0 or a%l]
        return max(logg(rts,t,l) for t in S)
    best=0
    for a in range(l):
        v=logg(rts,Fr(a),l)
        if v>0:  # lift greedily 3 levels
            cur=[a]; 
            for d in range(1,4):
                cand=[c+b*l**d for c in cur for b in range(l)] if l<400 else cur
                sc=[(logg(rts,Fr(c),l),c) for c in cand]; m=max(sc)[0]; v=max(v,m)
                cur=[c for s,c in sc if s==m][:2]
                if l>=400: break
        best=max(best,v)
    return best
def Cinf(y,eps):
    rts=KG.roots_cube(y,eps)
    t=sp.symbols('t')
    def poly(m):
        return sp.Poly(sp.expand(sp.prod([t-(sp.Rational(a.numerator,a.denominator)+sp.I*sp.Rational((-b if (m>>k)&1 else b).numerator,b.denominator)) for k,(a,b) in enumerate(rts)])),t,domain='QQ_I')
    pol={m:poly(m) for m in KG.MASKS}
    for a,b in itertools.combinations(KG.MASKS,2): assert (pol[a]-pol[b]).degree()<=4
    Gam=max(abs(complex(sp.N(((pol[a]-pol[b]).all_coeffs()[::-1]+[0]*5)[4]))) for a,b in itertools.combinations(KG.MASKS,2))
    P=badprimes(rts); lg=0
    for l in P:
        if l>200000: print("skip huge prime",l); continue
        lg+=localmax(rts,l)/2
    return Gam/math.exp(lg/2),Gam,P
if __name__=="__main__":
    yJ=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    for o in [(0,1,2,3),(0,2,3,1),(0,1,3,2)]:
        y=[yJ[i] for i in o]
        print("yJ order",o,Cinf(y,(1,-1,-1,1))[:2],flush=True)
    y2=[Fr(-25,6),Fr(2,3),Fr(18,5),Fr(-1,10)]
    print("y2",Cinf(y2,(1,-1,-1,1)),flush=True)
