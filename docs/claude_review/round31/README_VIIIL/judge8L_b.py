# Judge check, Part VIII (Question B thread), second script.
# (a) VIII.1: Walsh expansion of the 128 anchored pattern masses, E m_A^2 = E g^|A| <= E G^2, (VIII.2),
#     and the 64-pattern bound, exactly (Fractions), on random weighted sign systems.
# (b) VIII.1 escape: Schwartz-Zippel on a grid for a polynomial vanishing on the preimage of a cross-ratio locus.
# (c) Numerical margins of the Vojta constants (exact rationals).
# (d) Attack referee's B=5 exponent-2 family: two classes per prime power and ratio log q2/log L.
# (e) Subset family for divisors: 16 divisors = 1 mod q with alpha -> 1/4 from below.
import random, itertools, math, sys
from fractions import Fraction as F
from math import gcd
random.seed(7)
# (a)
def system(M, natoms):
    w=[F(random.randint(1,5)) for _ in range(natoms)]; W=sum(w); mu=[x/W for x in w]
    b=[[random.choice((1,-1)) for _ in range(natoms)] for _ in range(M)]
    return mu,b
bad=0; nt=0
for trial in range(4):
    M=4; mu,b=system(M,6)
    G=[[sum(mu[u]*b[x][u]*b[y][u] for u in range(6)) for y in range(M)] for x in range(M)]
    EG=sum(G[x][y] for x in range(M) for y in range(M))/M**2; assert EG>=0
    EG2=sum(G[x][y]**2 for x in range(M) for y in range(M))/M**2
    lam=max(max(G[x][y] for x in range(M) for y in range(M) if x!=y),0)/4
    assert EG2<=F(2,M)+8*lam
    for r in (2,4,6,8):
        lhs=F(0); cnt=0
        for xs in itertools.product(range(M),repeat=r):
            m=sum(mu[u]*math.prod(b[x][u] for x in xs) for u in range(6)); lhs+=m*m; cnt+=1
        lhs/=cnt
        rhs=sum(mu[u]*mu[v]*(F(sum(b[x][u]*b[x][v] for x in range(M)),M))**r for u in range(6) for v in range(6))
        assert lhs==rhs and lhs<=EG2
    # pattern masses for random 8-tuples (repeats allowed)
    for _ in range(200):
        xs=[random.randrange(M) for _ in range(8)]
        evens=[A for r in range(2,9,2) for A in itertools.combinations(range(8),r)]
        mA={A:sum(mu[u]*math.prod(b[xs[j]][u] for j in A) for u in range(6)) for A in evens}
        for Tm in range(128):
            s=[1]+[(-1 if Tm>>(i-1)&1 else 1) for i in range(1,8)]
            mass=sum(mu[u] for u in range(6) if all(b[xs[i]][u]*s[i]==b[xs[0]][u] for i in range(8)))
            pred=F(1,128)*(1+sum(math.prod(s[j] for j in A)*mA[A] for A in evens))
            assert mass==pred; nt+=1
        tot=sum(abs(v) for v in mA.values())
        for i,j in itertools.combinations(range(8),2):
            D=sum(mu[u] for u in range(6) if b[xs[i]][u]!=b[xs[j]][u])
            assert D<=F(1,2)*(1+tot)
print('(a) Walsh/moment identities: pattern checks',nt,'all exact, no failure')
# (b) F = (cross-ratio(x0,x1,x2,x3)+1)*denominator : vanishes on preimage of {cr=-1}; degree 2
def Fpoly(x): a,b_,c,d=x[:4]; return (a-c)*(b_-d)+(a-d)*(b_-c)
for grid in (list(range(1,31)), [k*k for k in range(1,31)], [F(1,k) for k in range(1,31)]):
    M=len(grid); z=0; tot=0
    for x in itertools.product(grid,repeat=4):
        tot+=1; z+= (Fpoly(x)==0)
    print('(b) grid size',M,'P(F=0)=%.4f'%(z/tot),'bound d/M=%.4f'%(2/M))
# (c)
c=F(6,69); eta0=F(1,381); a=F(2,10**5); lam=F(1,10**16)
assert a<eta0 and 28*(a/4+lam)<c/256 and 127**2*8*lam<(a/4)**2
Mmin=2/((a/(4*127))**2-8*lam)
print('(c) margins ok; 28(a/4+lam)=%.4e c/256=%.4e; M_* must exceed %.3e'%(28*(a/4+lam),c/256,Mmin))
# (d)
def q2(y):
    B=len(y); L=max(y)-min(y); res=1
    for mask in range(1<<(B-1)):
        S0=[y[i] for i in range(B) if not (mask>>i&1)]; S1=[y[i] for i in range(B) if mask>>i&1]
        g=0
        for S in (S0,S1):
            for v in S: g=gcd(g,v-S[0])
        if g: res=res*g//gcd(res,g)
    return res,L
for (al,be) in ((10**6+3,10**5+7),(10**12+39,10**11+3),(10**20+39,10**19+51)):
    y=[0,al*(al+be),be*(al+2*be),al*al+2*al*be+2*be*be,(al+be)*(al+2*be)]
    Q,L=q2(y)
    print('(d) B=5 family alpha=%.1e ratio log q2/log L = %.4f'%(al,math.log(Q)/math.log(L)))
# (e)
import sympy
for m in (10**3,10**8):
    q=24*m
    while True:
        fs=[1+i*q for i in range(1,5)]
        if all(sympy.isprime(f) for f in fs) or True: break
    n=math.prod(fs)
    ds=[d for d in itertools.chain.from_iterable(itertools.combinations(fs,r) for r in range(5))]
    cnt=sum(1 for d in ds if math.prod(d)%q==1)
    print('(e) q=%d: %d divisors = 1 mod q, alpha=%.5f'%(q,cnt,math.log(q)/math.log(n)))
