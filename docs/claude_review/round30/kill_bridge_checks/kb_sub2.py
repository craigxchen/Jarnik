# Fast sub-cluster constants: per prime pi, per sample t, valuations v_pi(z_j(t)) for all 8 points; subsets by min.
import itertools, math, sympy as sp
from fractions import Fraction as Fr
from kb_gauss import *
from kb_local import badprimes, vl
y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
rts=roots_cube(y); P=badprimes(rts)
def gprimes(l):
    if l==2: return [(1,1)]
    if l%4==3: return [(l,0)]
    a,b=sp.ntheory.factor_.primefactors(l) and None, None
    for x in range(1,int(l**0.5)+1):
        r=l-x*x; s=math.isqrt(r)
        if s*s==r: return [(x,s),(x,-s)]
def vpi_gint(z,pi):
    if z==(0,0): return 10**6
    n=pi[0]**2+pi[1]**2; v=0
    while True:
        x=z[0]*pi[0]+z[1]*pi[1]; yy=z[1]*pi[0]-z[0]*pi[1]
        if x%n or yy%n: return v
        z=(x//n,yy//n); v+=1
def vpi(val,pi,l):
    a,b=val; den=a.denominator*b.denominator//math.gcd(a.denominator,b.denominator)
    z=(int(a*den),int(b*den))
    e=2 if l==2 else (1 if l%4==1 else 1)
    vd=vl(Fr(den),l)*(2 if l==2 else 1)   # v_pi(den) for rational den
    return vpi_gint(z,pi)-vd
def samples(l):
    if l==2: return [Fr(a,2**j) for j in range(0,6) for a in range(2**(j+8)) if j==0 or a%2]
    if l==3: return [Fr(a,3**j) for j in range(0,4) for a in range(3**(j+5)) if j==0 or a%3]
    if l==5: return [Fr(a,5**j) for j in range(0,2) for a in range(5**(j+4)) if j==0 or a%5]
    return [Fr(a) for a in range(l)]
table={}
for l in P:
    pis=gprimes(l); rows=[]
    for t0 in samples(l):
        vals=values(rts,t0)
        rows.append([[vpi(v,pi,l) for v in vals] for pi in pis])
    table[l]=(pis,rows); print("done",l,len(rows),flush=True)
t=sp.symbols('t')
def poly(m):
    rr=[]
    for k,(a,b) in enumerate(rts):
        bb=-b if (m>>k)&1 else b
        rr.append(sp.Rational(a.numerator,a.denominator)+sp.I*sp.Rational(bb.numerator,bb.denominator))
    return sp.Poly(sp.expand(sp.prod([t-x for x in rr])),t,domain='QQ_I')
polys={m:poly(m) for m in MASKS}
def gam(a,b):
    co=(polys[a]-polys[b]).all_coeffs()[::-1]
    return abs(complex(sp.N(co[4]))) if len(co)>4 else 0.0
def Cinf(S):
    idx=[MASKS.index(m) for m in S]; lg=0; facs={}
    for l,(pis,rows) in table.items():
        best=max(sum(math.log(pi[0]**2+pi[1]**2)*min(r[k][j] for j in idx) for k,pi in enumerate(pis)) for r in rows)
        lg+=best/2
        if abs(best)>1e-9: facs[l]=math.exp(best/2)
    Gam=max(gam(a,b) for a,b in itertools.combinations(S,2))
    return Gam/math.exp(lg/2),facs,Gam
print("full", Cinf(MASKS))
best={}
for k in (5,6,7):
    for S in itertools.combinations(MASKS,k):
        c,f,G=Cinf(S)
        if k not in best or c<best[k][0]: best[k]=(c,[format(m,'08b')[::-1] for m in S],f)
for k in best: print(k,best[k])
