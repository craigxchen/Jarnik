# [J] independent exact check of the VI.4 four-point Pell families
import math, sys
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def N(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod(a,b):
    n=N(b); num=gmul(a,conj(b))
    q=((2*num[0]+n)//(2*n),(2*num[1]+n)//(2*n))
    r=(a[0]-gmul(q,b)[0],a[1]-gmul(q,b)[1]); return q,r
def ggcd(a,b):
    while b!=(0,0):
        _,r=gdivmod(a,b); a,b=b,r
    return a
def gexact(a,b):
    q,r=gdivmod(a,b); assert r==(0,0); return q
def Apow(v,k):
    x,y=v
    for _ in range(k): x,y=5*x+2*y,2*x+y
    return (x,y)
def run(seeds,rates,n):
    G=[Apow(v,r*n) for v,r in zip(seeds,rates)]
    Gb=[conj(g) for g in G]
    mi=(0,-1)
    pts=[gmul(gmul(Gb[0],Gb[1]),Gb[2]),
         gmul(mi,gmul(gmul(G[0],G[1]),Gb[2])),
         gmul(mi,gmul(gmul(G[0],Gb[1]),G[2])),
         gmul(mi,gmul(gmul(Gb[0],G[1]),G[2]))]
    raw=N(pts[0]); assert all(N(p)==raw for p in pts)
    g=pts[0]
    for p in pts[1:]: g=ggcd(g,p)
    q=[gexact(p,g) for p in pts]
    H=N(q[0]); assert all(N(p)==H for p in q)
    assert len(set(q))==4
    d2=max(N((a[0]-b[0],a[1]-b[1])) for a in q for b in q)
    # Delta = sqrt(d2/H); arc ~ sqrt(d2), R = sqrt(H)
    e=0.5*(math.log(d2)-math.log(H))/math.log(H)
    c13=math.exp(0.5*math.log(d2)-math.log(H)/6)
    c25=math.exp(0.5*math.log(d2)-math.log(H)/5)
    return N(g),H,e,c13,c25
seeds=[(3,1),(1,2),(2,-1)]
for rates in [(1,1,1),(6,7,7)]:
    print("rates",rates,"seeds",seeds)
    for n in ([5,10,15,20,25,30] if rates==(1,1,1) else [2,4,6,8,10,12]):
        Ng,H,e,c13,c25=run(seeds,rates,n)
        print(f" n={n:3d} N(gcd)={Ng} log10H={math.log10(H):.1f} exponent(log Delta/log H)={e:.5f} chord/R^(1/3)={c13:.5f} chord/R^(2/5)={c25:.5f}")
