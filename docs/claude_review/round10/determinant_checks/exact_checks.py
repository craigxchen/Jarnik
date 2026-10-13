# Exact checks (Gaussian-integer arithmetic) of the identities used in Thm 4.1 and Prop 5.5(i), and Lemma 5.4.
import itertools, random, math
from fractions import Fraction
import sympy as sp
random.seed(7)
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gpow(a,e):
    r=(1,0)
    for _ in range(e): r=gmul(r,a)
    return r
def gdivexact(a,b):  # a/b exact in Z[i]
    n=b[0]**2+b[1]**2; t=gmul(a,gconj(b)); assert t[0]%n==0 and t[1]%n==0; return (t[0]//n,t[1]//n)
def gdivisible(a,b):
    n=b[0]**2+b[1]**2; t=gmul(a,gconj(b)); return t[0]%n==0 and t[1]%n==0
def vpi(a,pi):
    if a==(0,0): return math.inf
    v=0
    while gdivisible(a,pi): a=gdivexact(a,pi); v+=1
    return v
def gaussian_prime(p):
    for x in range(1,p):
        y2=p-x*x
        if y2<=0: break
        y=int(math.isqrt(y2))
        if y*y==y2: return (x,y)
primes=[5,13,17,29]; N=1
for p in primes: N*=p
pts=[(x,y) for x in range(-int(math.isqrt(N)),int(math.isqrt(N))+1) for y in [int(math.isqrt(N-x*x))] if y*y==N-x*x and y>0]
pts=sorted(pts,key=lambda z: math.atan2(z[1],z[0]))
n=4
checks=0; facechecks=0; levelchecks=0
for trial in range(40):
    i0=random.randrange(len(pts)-n); Z=pts[i0:i0+n]
    T=2
    A=[k for k in itertools.product(range(-T,T+1),repeat=n-1)]
    A=[tuple(list(k)+[-sum(k)]) for k in A]
    A=random.sample(A, 40)
    m=3
    mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=n-1) if sum(a)==t]
    Mx=sp.Matrix([[sp.prod([k[j]**a[j] for j in range(n-1)]) for k in A] for a in mons])
    ns=Mx.nullspace()
    for vec in ns[:3]:
        den=sp.ilcm(*[sp.fraction(x)[1] for x in vec]); c=[int(x*den) for x in vec]
        if all(x==0 for x in c): continue
        D=max(sum(abs(t) for t in k) for k in A); D+=D%2
        def G(k):
            r=gpow((N,0),(D-sum(abs(t) for t in k))//2)
            for x,kx in enumerate(k):
                r=gmul(r, gpow(Z[x],kx) if kx>0 else gpow(gconj(Z[x]),-kx))
            return r
        F=(0,0)
        for ck,k in zip(c,A):
            g=G(k); F=(F[0]+ck*g[0],F[1]+ck*g[1])
        for p in primes:
            pi=gaussian_prime(p); e=1
            a=[vpi(z,pi) for z in Z]
            for k in A:  # valuation formula of (I1)
                assert vpi(G(k),pi)==e*D//2+sum(kk*aa for kk,aa in zip(k,a)); checks+=1
            supp=[(ck,k) for ck,k in zip(c,A) if ck!=0]
            vp=lambda x: (0 if x%p else 1+vp(x//p)) if x!=0 else math.inf
            levels=[vp(ck)+sum(kk*aa for kk,aa in zip(k,a)) for ck,k in supp]
            mu=min(levels)
            assert vpi(F,pi)>=e*D//2+mu
            # Prop 5.5(i): F/pi^{eD/2+mu} mod pi equals sigma^{eD/2}(N/p)^{D/2} sum_{Phi} c'_k sigma^{v_k} rho^k
            if F!=(0,0):
                lhs=F
                for _ in range(e*D//2+mu): lhs=gdivexact(lhs,pi)
                # residue field Z[i]/pi = F_p : map x+iy -> x + y*s where s^2=-1 mod p and pi | (s - i)... find s with pi | (s - i)
                s=[t for t in range(p) if gdivisible((t,-1),pi)][0]
                red=lambda g:(g[0]+g[1]*s)%p
                sigma=red(gconj(pi))
                y=[gdivexact(z, gmul(gpow(pi,aa),gpow(gconj(pi),e-aa))) for z,aa in zip(Z,a)]
                rho=[red(yy)*pow(sigma,-aa,p)%p for yy,aa in zip(y,a)]
                tot=0
                for (ck,k),lv in zip(supp,levels):
                    if lv==mu:
                        v=vp(ck); cp=ck//p**v
                        term=cp*pow(sigma,v,p)
                        for x,kx in enumerate(k): term=term*pow(rho[x],kx,p)
                        tot+=term
                rhs=pow(sigma,e*D//2,p)*pow((N//p)%p,D//2,p)*tot%p
                assert red(lhs)==rhs, (red(lhs),rhs); facechecks+=1
            # Lemma 5.4: level sums along <k,a_p> vanish when #levels <= m
            lv2={}
            for ck,k in supp: lv2.setdefault(sum(kk*aa for kk,aa in zip(k,a)),0); lv2[sum(kk*aa for kk,aa in zip(k,a))]+=ck
            if len(lv2)<=m:
                assert all(v==0 for v in lv2.values()); levelchecks+=1
print("valuation-formula checks:",checks," face-residue checks:",facechecks," level-sum checks:",levelchecks," all passed")
