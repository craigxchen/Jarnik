# Independent checks for astra_lemma sec 2.2, 2.3 (referee "baker").
import random, itertools, math, sys
from fractions import Fraction as F
from sympy import Matrix, isprime
from sympy.matrices.normalforms import hermite_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
import mpmath as mp

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gpow(a,n):
    r=(1,0)
    for _ in range(n): r=gmul(r,a)
    return r
def gprime(p):
    for a in range(1,int(p**.5)+1):
        b2=p-a*a; b=int(round(b2**.5))
        if b*b==b2 and b>0: return (a,b)
SPLIT=[p for p in range(5,2000) if isprime(p) and p%4==1]
PI={p:gprime(p) for p in SPLIT}
# Gaussian rationals as pairs of Fractions
def qmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def qinv(a):
    n=a[0]*a[0]+a[1]*a[1]; return (a[0]/n, -a[1]/n)
def qpow(a,n):
    if n<0: a=qinv(a); n=-n
    r=(F(1),F(0)); b=a
    while n:
        if n&1: r=qmul(r,b)
        b=qmul(b,b); n>>=1
    return r
def g_of(p):
    pi=PI[p]; return qmul((F(pi[0]),F(pi[1])), qinv((F(pi[0]),F(-pi[1]))))

def build(cols):
    """cols: dict p -> (e_p, [a_1..a_5]) ; returns z list"""
    z=[(1,0)]*5
    for p,(e,a) in cols.items():
        assert min(a)==0 and max(a)==e
        pi=PI[p]; pib=gconj(pi)
        for j in range(5):
            z[j]=gmul(z[j], gmul(gpow(pi,a[j]), gpow(pib,e-a[j])))
    return z

def analyse(cols, Smin=(), verbose=False):
    z=build(cols)
    assert len(set(z))==5
    H=z[0][0]**2+z[0][1]**2
    assert all(zz[0]**2+zz[1]**2==H for zz in z)
    W=math.log(H)
    assert abs(W-sum(e*math.log(p) for p,(e,a) in cols.items()))<1e-6
    zq=[(F(a),F(b)) for a,b in z]
    y=[qmul(zq[j],qinv(zq[4])) for j in range(5)]
    mp.mp.dps=int(W/2.3)+60
    yc=[mp.mpc(mp.mpf(t[0].numerator)/t[0].denominator, mp.mpf(t[1].numerator)/t[1].denominator) for t in y]
    Delta=max(abs(yc[i]-yc[j]) for i in range(5) for j in range(5))
    # columns
    S0=set(Smin)
    b={p:[a[j]-a[4] for j in range(4)] for p,(e,a) in cols.items()}
    # check y_j = eps_j prod g_p^{b_jp}, eps in mu4
    gp={p:g_of(p) for p in cols}
    for j in range(4):
        prod=(F(1),F(0))
        for p in cols: prod=qmul(prod,qpow(gp[p],b[p][j]))
        eps=qmul(y[j],qinv(prod))
        assert eps in [(F(1),F(0)),(F(-1),F(0)),(F(0),F(1)),(F(0),F(-1))], eps
    return z,y,yc,Delta,W,b,gp

def lattice_index(colvecs):
    M=Matrix(colvecs).T  # 4 x k
    if M.rank()<4: return None,None
    g=0
    for S in itertools.combinations(range(M.shape[1]),4):
        d=M[:,list(S)].det(); g=math.gcd(g,int(d))
    Hn=hermite_normal_form(M)
    I2=abs(int(Hn.det())) if Hn.shape==(4,4) else None
    assert I2==g, (I2,g)
    return g,Hn

def short_dual(Hn,I):
    P=Matrix(Hn)  # columns basis of L
    D=(P.T).inv()*I  # columns: I * basis of L*
    D=D.applyfunc(lambda t:int(t))
    dm=DomainMatrix([[ZZ(int(D[i,j])) for i in range(4)] for j in range(4)],(4,4),ZZ)  # rows = basis vectors
    red=dm.lll().to_Matrix()
    rows=[list(red.row(i)) for i in range(4)]
    best=None
    R=3
    for co in itertools.product(range(-R,R+1),repeat=4):
        if not any(co): continue
        w=[sum(co[i]*rows[i][k] for i in range(4)) for k in range(4)]
        n=max(abs(t) for t in w)
        if best is None or n<best[0]: best=(n,w)
    n,w=best
    v=[F(int(t),I) for t in w]
    return v

def gauss_ratio_norm_check(u_exps,gp):
    pass

def main_check(cols, tag):
    z,y,yc,Delta,W,b,gp=analyse(cols)
    I,Hn=lattice_index([b[p] for p in cols])
    if I is None: return None
    assert I <= 24*(W/math.log(5))**4
    v=short_dual(Hn,I)
    vinf=max(abs(t) for t in v)
    assert vinf <= F(2)*1 and float(vinf) <= I**(-0.25)+1e-12, (vinf,I)   # sharper Minkowski t=I^{-1/4}
    # v in L*
    k={p:sum(v[j]*b[p][j] for j in range(4)) for p in cols}
    assert all(t.denominator==1 for t in k.values())
    assert any(t!=0 for t in k.values())
    den=1
    for t in v: den=den*t.denominator//math.gcd(den,t.denominator)
    m=den
    assert I % m ==0
    c=[int(m*t) for t in v]
    c1=sum(abs(t) for t in c); assert c1<=8*I
    twoh=sum(abs(int(k[p]))*math.log(p) for p in cols)
    assert twoh <= 4*float(vinf)*W+1e-9
    # u exact
    u=(F(1),F(0))
    for p in cols: u=qmul(u,qpow(gp[p],int(k[p])))
    yc_c=(F(1),F(0))
    for j in range(4): yc_c=qmul(yc_c,qpow(y[j],c[j]))
    lhs=qpow(yc_c,4); rhs=qpow(u,4*m)
    assert lhs==rhs
    uc=mp.mpc(mp.mpf(u[0].numerator)/u[0].denominator, mp.mpf(u[1].numerator)/u[1].denominator)
    val=abs(uc**(4*m)-1)
    assert val>0
    assert val <= 4*c1*Delta*(1+mp.mpf(10)**(-30))
    hu=twoh/2
    ratio=float(-mp.log(val))/((1+hu)*math.log(8*m))
    lhs2=float(-mp.log(Delta))/W
    return dict(tag=tag,I=I,m=m,vinf=float(vinf),W=W,h=hu,ratio=ratio,expo=lhs2, I14=I**-0.25*math.log(2*I))

def rand_cols(nprimes, E, mult=1, rng=random):
    ps=rng.sample(SPLIT[:60],nprimes)
    cols={}
    for p in ps:
        e=rng.randint(1,E)
        while True:
            a=[rng.randint(0,e) for _ in range(5)]
            if min(a)==0 and max(a)==e: break
        cols[p]=(e*mult,[t*mult for t in a])
    return cols

if __name__=="__main__":
    rng=random.Random(31)
    res=[]
    for trial in range(300):
        np_=rng.randint(4,7); E=rng.randint(1,4); mult=rng.choice([1,1,1,2,3,5])
        cols=rand_cols(np_,E,mult,rng)
        if len(set(build(cols)))<5: continue
        r=main_check(cols,'rand')
        if r: res.append(r)
    print("random full-rank tuples checked:",len(res))
    print("max ratio -log|u^4m-1|/((1+h)log 8m):",max(r['ratio'] for r in res))
    print("max I:",max(r['I'] for r in res))
    # structured: y_j = x_j^N  (all columns multiplied by N)
    base=None
    for trial in range(50):
        cols=rand_cols(4,2,1,rng)
        z=build(cols)
        if len(set(z))<5: continue
        if lattice_index([[a[j]-a[4] for j in range(4)] for p,(e,a) in cols.items()])[0] is None: continue
        base=cols;break
    for N in [1,2,3,5,8,13,21]:
        cols={p:(e*N,[t*N for t in a]) for p,(e,a) in base.items()}
        r=main_check(cols,'pow%d'%N)
        print(r)
