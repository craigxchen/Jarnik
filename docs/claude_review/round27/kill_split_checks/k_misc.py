# Lemma 3.5 exact Vandermonde identity, Theorem 3.6 cross product, phase-0 k-tuples. Own code.
import itertools, random, math, sympy as sp
from math import gcd
from functools import reduce
random.seed(7)
Bs=[frozenset(s) for r in range(1,5) for s in itertools.combinations(range(5),r)]
def nB_of(u):
    # n_B = product over layers (p,s) of type B of p
    n={B:1 for B in Bs}
    ps=set()
    for x in u: ps|=set(sp.factorint(x))
    for p in ps:
        v=[sp.multiplicity(p,x) for x in u]
        for s in range(1,max(v)+1):
            A=frozenset(i for i in range(5) if v[i]>=s)
            if 0<len(A)<5: n[A]*=p
    return n
ok36=0; okV=0; N=200
for trial in range(N):
    while True:
        u=[random.randint(2,10**6) for _ in range(5)]
        g=reduce(gcd,u); u=[x//g for x in u]
        if len(set(u))==5: break
    n=nB_of(u)
    assert all(math.prod(n[B] for B in Bs if a in B)==u[a] for a in range(5))
    # Thm 3.6 for every triple
    good=True
    for tau in itertools.combinations(range(5),3):
        de=[i for i in range(5) if i not in tau]
        def N_(T): return math.prod(n[frozenset(T)|frozenset(R)] for r in range(3) for R in itertools.combinations(de,r) if 0<len(frozenset(T)|frozenset(R))<5)
        a,b,c=tau
        x=[N_([a]),N_([b]),N_([c])]; y=[N_([b,c]),N_([c,a]),N_([a,b])]
        cr=(x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])
        gab=lambda i,j: math.prod(n[B] for B in Bs if i in B and j in B)
        t=lambda i,j: (u[i]-u[j])//gab(i,j)
        if cr!=(t(b,c),t(c,a),t(a,b)): good=False
    ok36+=good
    # Vandermonde: V(y)=prod(u_a-u_b)/prod u_a^2 == prod t_ab * prod n_B^{-|B|(5-|B|)/2}
    lhs=sp.Rational(math.prod(u[i]-u[j] for i,j in itertools.combinations(range(5),2)), math.prod(x*x for x in u))
    rhs=sp.Rational(math.prod((u[i]-u[j])//math.prod(n[B] for B in Bs if i in B and j in B) for i,j in itertools.combinations(range(5),2)))
    rhs*=sp.Rational(1, math.prod(n[B]**(len(B)*(5-len(B))//2) for B in Bs))
    okV+=(lhs==rhs)
print("Thm 3.6 cross-product identity holds on",ok36,"/",N,"random primitive 5-tuples")
print("Vandermonde exact formula V(y)=prod t * prod n_B^{-f(|B|)} holds on",okV,"/",N)
# sum_B x_B f(|B|) = 40 on cut-balanced profiles (random)
for _ in range(3):
    x={}
    for B in Bs:
        if 0 in B: x[B]=sp.Rational(random.randint(0,100),100); x[frozenset(range(5))-B]=1-x[B]
    print("sum x_B |B|(5-|B|)/2 =",sum(x[B]*len(B)*(5-len(B)) for B in Bs)/2)
# phase 0 for m-czech: n=(x^2-1)(x^2-4)(x^2-16)(x^2-49), divisors F_0..F_4 at Pell points z^2=56x^2+65
def F(xx,zz):
    a,b,c=1,2,4; d=7; S=1+4+16+49; m=56
    return [(xx-a)*(xx-b)*(xx-c)*(xx+d), xx**4-S*xx**2//2+m+2*xx*zz,(xx**2-b*b)*(xx**2-c*c),(xx**2-a*a)*(xx**2-c*c),(xx**2-a*a)*(xx**2-b*b)]
pts=[]
X,Z=None,None
for xx in range(1,10**6):
    zz2=56*xx*xx+65; zz=math.isqrt(zz2)
    if zz*zz==zz2: pts.append((xx,zz))
for (xx,zz) in pts:
    if xx<8: continue
    n=(xx*xx-1)*(xx*xx-4)*(xx*xx-16)*(xx*xx-49)
    D=F(xx,zz)
    if not all(dd>0 and n%dd==0 for dd in D): print("  x",xx,"not all divisors"); continue
    s=[dd+n//dd for dd in D]; k=[si-min(s) for si in s]
    z=[dd-n//dd for dd in D]
    idt=all(z[j]**2-z[0]**2==2*(s[j]-s[0])*s[0]+(s[j]-s[0])**2 for j in range(5))
    print("  m-czech x=",xx,"k (sorted)=",sorted(k),"identity:",idt, "C=",(max(D)-math.isqrt(n))/n**0.25)
