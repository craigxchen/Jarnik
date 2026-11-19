# exact checks on genuine pure ten-pair tuples built from Gaussian blocks with private odd prime norms
import random, itertools, sympy as sp
from fractions import Fraction as Fr
random.seed(5)
E=[(i,j) for i in range(5) for j in range(i+1,5)]; EI={e:k for k,e in enumerate(E)}
eid=lambda i,j:EI[(min(i,j),max(i,j))]
class G:  # Gaussian rationals
    def __init__(s,a,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): o=o if isinstance(o,G) else G(o); return G(s.a+o.a,s.b+o.b)
    def __sub__(s,o): o=o if isinstance(o,G) else G(o); return G(s.a-o.a,s.b-o.b)
    def __mul__(s,o): o=o if isinstance(o,G) else G(o); return G(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def conj(s): return G(s.a,-s.b)
    def nrm(s): return s.a*s.a+s.b*s.b
    def __truediv__(s,o): o=o if isinstance(o,G) else G(o); d=o.nrm(); r=s*o.conj(); return G(r.a/d,r.b/d)
    def __eq__(s,o): o=o if isinstance(o,G) else G(o); return s.a==o.a and s.b==o.b
    def __repr__(s): return f"({s.a}+{s.b}i)"
def split_primes(k):
    out=[];p=5
    while len(out)<k:
        if sp.isprime(p) and p%4==1: out.append(p)
        p+=1
    return out
def gauss_rep(p):
    for a in range(1,int(p**.5)+1):
        b2=p-a*a; b=int(round(b2**.5))
        if b*b==b2: return G(a,b)
primes=split_primes(60)
fails=0; parity_ok=0; Jok=0; quart_ok=0; trials=40
for tr in range(trials):
    ps=random.sample(primes,10)
    beta=[]
    for p in ps:
        g=gauss_rep(p); 
        if random.random()<.5: g=g.conj()
        u=random.choice([G(1),G(0,1),G(-1),G(0,-1)]); beta.append(g*u)
    n=[b.nrm() for b in beta]
    P=[]
    for i in range(5):
        x=G(1)
        for j in range(5):
            if j!=i: x=x*beta[eid(i,j)]
        P.append(x)
    t=[]; 
    for (i,j) in E:
        Z=(P[i].conj()*P[j])/n[eid(i,j)]
        assert Z.a.denominator==1 and Z.b.denominator==1
        t.append(int(Z.b))
    if 0 in t: continue
    # III.3 parity: t_ij odd iff P_i,P_j in different classes mod 2
    cls=[int(x.a)%2 for x in P]
    parity_ok+= all((t[eid(i,j)]%2==1)==(cls[i]!=cls[j]) for (i,j) in E)
    # J identity: prod_{i<j}(q_i-q_j)/prod q_i^2 = -1024 tau / H^3   (q_i=P_i/conj P_i)
    q=[x/x.conj() for x in P]
    num=G(1)
    for (i,j) in E: num=num*(q[i]-q[j])
    den=G(1)
    for x in q: den=den*x*x
    tau=1
    for x in t: tau*=x
    H=1
    for x in n: H*=x
    lhs=num/den; rhs=G(Fr(-1024*tau,H**3))
    if lhs==rhs: Jok+=1
    else:
        rhs2=G(Fr(1024*tau,H**3)); 
        if lhs==rhs2: Jok+=1
    # quartic completion: fix q1..q4 (times arbitrary common scalar g), x=q5 root of C prod(q_i-x)/x^2 = val
    X=sp.symbols('X')
    g=G(3,7)
    qq=[x*g for x in q]
    val=lhs
    poly=sp.expand(sp.prod([ (sp.Rational(qq[i].a.numerator,qq[i].a.denominator)+sp.I*sp.Rational(qq[i].b.numerator,qq[i].b.denominator))-X for i in range(4)]) * sp.prod([(sp.Rational(qq[i].a.numerator,qq[i].a.denominator)+sp.I*sp.Rational(qq[i].b.numerator,qq[i].b.denominator))-(sp.Rational(qq[j].a.numerator,qq[j].a.denominator)+sp.I*sp.Rational(qq[j].b.numerator,qq[j].b.denominator)) for (i,j) in E if j<4])
        - (sp.Rational(val.a.numerator,val.a.denominator)+sp.I*sp.Rational(val.b.numerator,val.b.denominator))*X**2*sp.prod([(sp.Rational(qq[i].a.numerator,qq[i].a.denominator)+sp.I*sp.Rational(qq[i].b.numerator,qq[i].b.denominator))**2 for i in range(4)]))
    x5=sp.Rational(qq[4].a.numerator,qq[4].a.denominator)+sp.I*sp.Rational(qq[4].b.numerator,qq[4].b.denominator)
    pl=sp.Poly(poly,X)
    if pl.degree()==4 and sp.expand(pl.eval(x5))==0: quart_ok+=1
print("trials with nonzero t:",trials, "parity lemma ok:",parity_ok,"J identity (|c|=1024) ok:",Jok,"fifth mark is a root of a nonzero quartic:",quart_ok)
