# independent library: Gaussian points, arcs, twin classes (adversarial check of lane clothes)
import math, sympy
from math import gcd, isqrt
def gprime(p):  # pi with |pi|^2=p, p=1 mod 4
    r=sympy.sqrt_mod(-1,p); a,b=p,r
    while b*b>p: a,b=b,a%b
    c=isqrt(p-b*b); assert b*b+c*c==p; return (b,c)
def gmul(u,v): return (u[0]*v[0]-u[1]*v[1], u[0]*v[1]+u[1]*v[0])
def points(N):
    f=sympy.factorint(N); cur=[(1,0)]
    for p,e in f.items():
        if p==2:
            for _ in range(e): cur=[gmul(z,(1,1)) for z in cur]
        elif p%4==3:
            assert e%2==0; cur=[(z[0]*p**(e//2),z[1]*p**(e//2)) for z in cur]
        else:
            pi=gprime(p); pb=(pi[0],-pi[1]); new=[]
            for z in cur:
                for a in range(e+1):
                    w=z
                    for _ in range(a): w=gmul(w,pi)
                    for _ in range(e-a): w=gmul(w,pb)
                    new.append(w)
            cur=new
    S=set()
    for z in cur:
        for u in [(1,0),(0,1),(-1,0),(0,-1)]: S.add(gmul(z,u))
    S=sorted(S); assert all(x*x+y*y==N for x,y in S); return S
def arcmax(N,C):
    P=points(N); R=math.sqrt(N); th=sorted(math.atan2(y,x) for x,y in P)
    n=len(th); th2=th+[t+2*math.pi for t in th]; best=1; j=0
    lim=C*math.sqrt(R)/R
    for i in range(n):
        if j<i: j=i
        while j+1<i+n and th2[j+1]-th2[i]<=lim+1e-12: j+=1
        best=max(best,j-i+1)
    return best
def twinmax(N,lam,P=None):
    P=P or points(N); Q=lam*N**0.25; best=(1,None,None)
    for z1 in P:
        cnt={}
        for z in P:
            if z==z1: continue
            G=gcd(z[0]-z1[0],z[1]-z1[1])
            for d in sympy.divisors(G):
                d=int(d)
                if d>=Q and d&1 and gcd(d,N)==1: cnt[d]=cnt.get(d,0)+1
        for d,v in cnt.items():
            if v+1>best[0]: best=(v+1,d,z1)
    return best
def classes(P,q):
    D={}
    for z in P: D.setdefault((z[0]%q,z[1]%q),[]).append(z)
    return D
