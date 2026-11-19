# Index of all 5-subsets on the record's explicit genuine clusters (growth/families.md), plus hypothesis checks.
import math, sys
from itertools import combinations
from sympy import factorint
from lat import anchored_index
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gdivexact(a,b):
    n=b[0]**2+b[1]**2
    x=a[0]*b[0]+a[1]*b[1]; y=a[1]*b[0]-a[0]*b[1]
    if x%n or y%n: return None
    return (x//n,y//n)
def gprime(p):
    for x in range(1,int(math.isqrt(p))+1):
        y2=p-x*x; y=math.isqrt(y2)
        if y*y==y2: return (x,y)
def ggcd(a,b):
    while b!=(0,0):
        n=b[0]**2+b[1]**2
        x=a[0]*b[0]+a[1]*b[1]; y=a[1]*b[0]-a[0]*b[1]
        qx=round(x/n); qy=round(y/n)
        r=(a[0]-(qx*b[0]-qy*b[1]), a[1]-(qx*b[1]+qy*b[0]))
        a,b=b,r
    return a
def points(n):
    pts=[]
    r=math.isqrt(n)
    for x in range(-r,r+1):
        y2=n-x*x; y=math.isqrt(y2)
        if y*y==y2:
            pts.append((x,y))
            if y: pts.append((x,-y))
    return pts
def best_cluster(n,k):
    pts=sorted(points(n),key=lambda p: math.atan2(p[1],p[0]))
    L=len(pts); R=math.sqrt(n); best=None
    for i in range(L):
        j=(i+k-1)%L
        a=(math.atan2(pts[j][1],pts[j][0])-math.atan2(pts[i][1],pts[i][0]))%(2*math.pi)
        C=a*R/math.sqrt(R)
        if best is None or C<best[0]: best=(C,[pts[(i+t)%L] for t in range(k)])
    return best
def analyse(name,pts):
    g=pts[0]
    for p in pts[1:]: g=ggcd(g,p)
    red=[gdivexact(p,g) for p in pts]
    n=red[0][0]**2+red[0][1]**2
    R=math.sqrt(n)
    f=factorint(n)
    assert 2 not in f, f
    primes=sorted(f)
    rows=[]
    for z in red:
        row=[]
        for p in primes:
            pi=gprime(p); a=0; w=z
            while True:
                t=gdivexact(w,pi)
                if t is None: break
                w=t; a+=1
            row.append(a)
        rows.append(row)
    # unit class: z / prod pi^a pibar^(e-a) is a unit
    units=[]
    for z,row in zip(red,rows):
        P=(1,0)
        for p,a in zip(primes,row):
            pi=gprime(p); pib=(pi[0],-pi[1])
            for _ in range(a): P=gmul(P,pi)
            for _ in range(f[p]-a): P=gmul(P,pib)
        units.append(gdivexact(z,P))
    W=math.log(n)
    ang=[math.atan2(z[1],z[0]) for z in red]
    m=len(red)
    ang=[(a-ang[0]+math.pi)%(2*math.pi)-math.pi for a in ang]; width=max(ang)-min(ang)
    seps=[sum(abs(rows[i][t]-rows[j][t])*math.log(primes[t]) for t in range(len(primes))) for i,j in combinations(range(m),2)]
    idx={}
    for T in combinations(range(m),5):
        I=anchored_index(rows,T); idx[I]=idx.get(I,0)+1
    print(f"{name}: m={m}, reduced N={n}={f}, W={W:.2f}, arc C={width*R/math.sqrt(R):.4f}, angular width/(e^(-W/4)/2)={width/(0.5*math.exp(-W/4)):.2f}")
    print(f"   unit classes {sorted(set(units))}, min pair sep/(W/2)={min(seps)/(W/2):.3f}, #pairs with sep<=W/2: {sum(s<=W/2 for s in seps)}/{len(seps)}")
    print(f"   affine rank of rows = {affrank(rows)} (m-1={m-1}); 5-subset index distribution {dict(sorted(idx.items()))}")
def affrank(rows):
    from sympy import Matrix
    return Matrix([[a-b for a,b in zip(r,rows[0])] for r in rows[1:]]).rank()
for n,k,name in [(32988780417125,8,"8-pt record C=5.674"),(235370525,7,"7-pt record C=5.27"),(537606725,8,"golden n=4"),(276108509,3,"3-pt"),(412333805,4,"4-pt")]:
    C,pts=best_cluster(n,k)
    if k>=5: analyse(name,pts)
    else: print(name,"C=",C)
