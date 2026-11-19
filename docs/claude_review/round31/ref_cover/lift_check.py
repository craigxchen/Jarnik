# (b),(c): toric height = W computed directly from place-by-place valuations; lift x with y=eps f_P(x) exactly.
import random, sys
from itertools import combinations
from math import gcd, log
from fractions import Fraction
import sympy as sp
sys.path.insert(0,'.')
from gauss import *

def qi_div(a,b):  # Gaussian rational as pair of Fractions
    n=gnorm(b); num=gmul(a,gconj(b)); return (Fraction(num[0],n),Fraction(num[1],n))
def qmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def qinv(a):
    n=a[0]**2+a[1]**2; return (a[0]/n, -a[1]/n)
def qpow(a,k):
    if k<0: a=qinv(a); k=-k
    r=(Fraction(1),Fraction(0))
    for _ in range(k): r=qmul(r,a)
    return r
def qval(q,pi):
    # valuation at pi of Gaussian rational q: write q = num/den with den rational integer
    from math import lcm
    d=lcm(q[0].denominator,q[1].denominator)
    num=(int(q[0]*d),int(q[1]*d))
    return val(num,pi)-val((d,0),pi)

random.seed(1)
trials=0; nonfull=0
for t in range(300):
    ps=random.sample(SPLIT[:12], random.randint(2,6))
    es={p:random.randint(1,3) for p in ps}
    # random allocations with min 0, max e_p
    while True:
        A={}
        for p in ps:
            col=[random.randint(0,es[p]) for _ in range(5)]
            i,j=random.sample(range(5),2); col[i]=0; col[j]=es[p]
            A[p]=col
        us=[random.choice(UNITS) for _ in range(5)]
        z=[]
        for i in range(5):
            w=us[i]
            for p in ps:
                w=gmul(w,gmul(gpow(PI[p],A[p][i]),gpow(gconj(PI[p]),es[p]-A[p][i])))
            z.append(w)
        if len(set(z))==5: break
    H=1
    for p in ps: H*=p**es[p]
    assert all(gnorm(w)==H for w in z)
    g=z[0]
    for w in z[1:]: g=ggcd(g,w)
    assert gnorm(g)==1, "not primitive"
    y=[qi_div(z[i],z[4]) for i in range(5)]
    # toric height: prod over places of Q(i) of max_{i,j}|y_i/y_j|_w ^{1/2} (each split place local degree 1, global degree 2)
    logHB=0.0; HB=Fraction(1)
    for p in ps:
        for pi in (PI[p],gconj(PI[p])):
            vs=[qval(yy,pi) for yy in y]
            spread=max(vs)-min(vs)       # max |y_i/y_j|_pi = p^{spread}
            HB*=Fraction(p)**spread
    # archimedean: all |y_i|=1 -> factor 1. HB is prod p^{spread}; absolute height = HB^{1/2}
    assert HB==Fraction(H)**2, (HB,H)   # i.e. absolute toric height = H exactly
    # columns
    cols=[[A[p][j]-A[p][4] for j in range(4)] for p in ps]
    M=sp.Matrix(cols).T
    if M.rank()<4: nonfull+=1; continue
    minors=[abs(M[:,list(c)].det()) for c in combinations(range(len(ps)),4)]
    I=0
    for m in minors: I=gcd(I,int(m))
    # HNF basis of column lattice
    from sympy.matrices.normalforms import hermite_normal_form
    Hn=hermite_normal_form(M)
    P=Hn[:, :4] if Hn.shape[1]>=4 else Hn
    assert abs(P.det())==I, (P, I)
    Pinv=P.inv()
    x=[]
    hs={}
    for k,p in enumerate(ps):
        h=Pinv*sp.Matrix(cols[k]); assert all(v.is_integer for v in h); hs[p]=[int(v) for v in h]
    gp={p:qi_div(PI[p],gconj(PI[p])) for p in ps}
    for kk in range(4):
        w=(Fraction(1),Fraction(0))
        for p in ps: w=qmul(w,qpow(gp[p],hs[p][kk]))
        x.append(w)
    for j in range(4):
        fx=(Fraction(1),Fraction(0))
        for kk in range(4): fx=qmul(fx,qpow(x[kk],int(P[j,kk])))
        eps=qmul(y[j],qinv(fx))
        assert eps in [(Fraction(a),Fraction(b)) for a,b in UNITS], eps
    trials+=1
print("full-rank tuples checked:",trials,"rank<4:",nonfull)
