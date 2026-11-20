# Exact check of the point-residue design rho_i = s0*omega_i/conj(omega_i) mod b, b odd prime power <= 2m,
# against the record's D_u(c) (round5/construct.md Def 1.2): D_u | A - u*conj(A), (1+i)D_u | A-u conj A,
# D_u <= sqrt2 L0^{s/2}; pairs: D_1 | Im(omega_i conj omega_j), |sin(beta_i-beta_j)|>=1/L0.
import math, random, itertools
from sympy import isprime, primerange
def gmul(a,b,M=None):
    r=(a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    return (r[0]%M,r[1]%M) if M else r
def gconj(a): return (a[0],-a[1])
def gpow(a,e,M):
    r=(1,0); a=(a[0]%M,a[1]%M)
    while e:
        if e&1: r=gmul(r,a,M)
        a=gmul(a,a,M); e>>=1
    return r
def ginv(a,M):  # a unit mod M (norm coprime to M)
    n=(a[0]*a[0]+a[1]*a[1])%M; ni=pow(n,-1,M); c=gconj(a); return (c[0]*ni%M,c[1]*ni%M)
def gauss_prime(l):
    for x in range(1,int(math.isqrt(l))+1):
        y2=l-x*x; y=math.isqrt(y2)
        if y*y==y2:
            a,b=max(x,y),min(x,y)   # argument in (0,pi/4)
            return (a,b)
rng=random.Random(3)
for m in [20,40,60]:
    ells=[]; p=2*m+1
    while len(ells)<m:
        if isprime(p) and p%4==1: ells.append(p)
        p+=1
    L0=max(ells); om=[gauss_prime(l) for l in ells]; beta=[math.atan2(b,a) for a,b in om]
    assert all(0<x<math.pi/4 for x in beta)
    moduli=[]  # odd prime powers <= 2m grouped by prime
    for r in primerange(3,2*m+1):
        pw=[]; b=r
        while b<=2*m: pw.append(b); b*=r
        moduli.append(pw)
    s0=gmul(gauss_prime(1000033),gauss_prime(1000037))  # norm a product of primes > 2m, as for N
    def rho(i,b): return gmul(gmul(s0,om[i],b), ginv(gconj(om[i]),b), b)
    units={1:(1,0),-1:(-1,0),'i':(0,1),'-i':(0,-1)}
    def D(c,u):
        Dv=1
        for pw in moduli:
            best=1
            for b in pw:
                prod=(1,0)
                for i,ci in c.items():
                    r=rho(i,b); prod=gmul(prod, gpow(r,ci,b) if ci>0 else gpow(ginv(r,b),-ci,b), b)
                if prod==((u[0])%b,(u[1])%b): best=b
            Dv*=best
        return Dv
    def Aof(c):
        A=(1,0)
        for i,ci in c.items():
            f=om[i] if ci>0 else gconj(om[i])
            for _ in range(abs(ci)): A=gmul(A,f)
        return A
    ok=True; minsin=1; maxratio=0; npair=0
    for i,j in itertools.combinations(range(m),2):
        c={i:1,j:-1}; A=Aof(c); npair+=1
        im=gmul(om[i],gconj(om[j]))[1]
        minsin=min(minsin,abs(math.sin(beta[i]-beta[j]))*L0)
        for name,u in units.items():
            d=D(c,u); uA=gmul(u,gconj(A)); diff=(A[0]-uA[0],A[1]-uA[1])
            ok&= diff!=(0,0) and diff[0]%d==0 and diff[1]%d==0 and (diff[0]//d+diff[1]//d)%2==0
            if name==1: ok&= im%d==0 and d<=abs(im)<=L0*abs(beta[i]-beta[j])+1e-9
    nlam=0
    for _ in range(400):
        s=rng.choice([4,6,8]); idx=rng.sample(range(m),s)
        c={}
        for t,i in enumerate(idx): c[i]=c.get(i,0)+(1 if t<s//2 else -1)
        A=Aof(c)
        for u in units.values():
            d=D(c,u); uA=gmul(u,gconj(A)); diff=(A[0]-uA[0],A[1]-uA[1])
            ok&= diff!=(0,0) and diff[0]%d==0 and diff[1]%d==0 and (diff[0]//d+diff[1]//d)%2==0 and d<=math.sqrt(2)*L0**(s/2)
            maxratio=max(maxratio,d/(math.sqrt(2)*L0**(s/2)))
        nlam+=1
    print("m=%d L0=%d pairs=%d (x4 units) chars=%d (x4): all checks %s; min L0|sin(beta_i-beta_j)| = %.3f (>=1); max D_u/(sqrt2 L0^(s/2)) = %.2e"%(m,L0,npair,nlam,ok,minsin,maxratio))
