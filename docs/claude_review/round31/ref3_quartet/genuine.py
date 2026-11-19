"""Referee check of (III.1) on genuine lattice-point clusters from growth/families.md and the scan data.
Independent code: Gaussian factorisation, exact integers only (angles in float only to select windows)."""
import sys, math
from itertools import combinations, product
from math import gcd, isqrt, prod, lcm
from sympy import factorint

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gpow(a,k):
    r=(1,0)
    for _ in range(k): r=gmul(r,a)
    return r
def gdivexact(a,b):
    n=gnorm(b); c=gmul(a,gconj(b))
    if c[0]%n or c[1]%n: return None
    return (c[0]//n,c[1]//n)
def two_sq(p):
    for a in range(1,isqrt(p)+1):
        b2=p-a*a; b=isqrt(b2)
        if b*b==b2: return (a,b)
def vpi(z,pi):
    k=0
    while True:
        q=gdivexact(z,pi)
        if q is None: return k
        z=q; k+=1
UNITS=[(1,0),(0,1),(-1,0),(0,-1)]

def points(n):
    f=factorint(n); assert all(p%4==1 for p in f), f
    pis={p:two_sq(p) for p in f}
    pts=[]
    for ex in product(*[range(e+1) for e in f.values()]):
        z=(1,0)
        for (p,e),a in zip(f.items(),ex):
            z=gmul(z,gmul(gpow(pis[p],a),gpow(gconj(pis[p]),e-a)))
        for u in UNITS: pts.append(gmul(u,z))
    assert all(gnorm(z)==n for z in pts)
    return sorted(set(pts),key=lambda z: math.atan2(z[1],z[0]))

def best_window(pts,k,n):
    ang=[math.atan2(z[1],z[0]) for z in pts]; L=len(pts); best=None
    for s in range(L):
        e=(s+k-1)%L; span=(ang[e]-ang[s])%(2*math.pi)
        if best is None or span<best[0]: best=(span,s)
    span,s=best
    return [pts[(s+j)%L] for j in range(k)], span*n**0.25

def ggcd(a,b):
    while b!=(0,0):
        # Euclid in Z[i]
        n=gnorm(b); c=gmul(a,gconj(b))
        q=(round(c[0]/n) if abs(c[0])<2**52 else (2*c[0]+n)//(2*n), (2*c[1]+n)//(2*n))
        q=((2*c[0]+n)//(2*n),(2*c[1]+n)//(2*n))
        r=(a[0]-gmul(q,b)[0],a[1]-gmul(q,b)[1]); a,b=b,r
    return a

def analyse(cluster):
    g=cluster[0]
    for z in cluster[1:]: g=ggcd(g,z)
    zs=[gdivexact(z,g) for z in cluster]
    N=gnorm(zs[0]); f=factorint(N); assert all(p%4==1 for p in f)
    pis={p:two_sq(p) for p in f}
    data=[]
    for z in zs:
        a={p:vpi(z,pis[p]) for p in f}
        assert all(vpi(z,gconj(pis[p]))==f[p]-a[p] for p in f)
        base=(1,0)
        for p in f: base=gmul(base,gmul(gpow(pis[p],a[p]),gpow(gconj(pis[p]),f[p]-a[p])))
        eps=next(u for u in UNITS if gmul(u,base)==z)
        P=(1,0)
        for p in f: P=gmul(P,gpow(pis[p],a[p]))
        data.append((z,a,eps,P))
    return N,f,data

def det(u,v): return u[0]*v[1]-u[1]*v[0]

def check(name,n,k):
    pts=points(n); cl,C=best_window(pts,k,n)
    N,f,data=analyse(cl)
    W=math.log(N)
    nq=nsame=0; bad=0; maxratio=0
    for q in combinations(range(k),4):
        Ps=[data[i][3] for i in q]; ex=[data[i][1] for i in q]
        g={}; t={}
        for i,j in combinations(range(4),2):
            g[i,j]=prod(p**min(ex[i][p],ex[j][p]) for p in f)
            D=det(Ps[i],Ps[j]); assert D%g[i,j]==0 and D!=0
            t[i,j]=D//g[i,j]
        raw=[g[0,1]*g[2,3],g[0,2]*g[1,3],g[0,3]*g[1,2]]
        G=gcd(*raw); A,B,Cc=[r//G for r in raw]
        assert gcd(A,B)==gcd(B,Cc)==gcd(A,Cc)==1 and N%(A*B*Cc)==0
        # matching rule h3-h2 on the low/high pairing
        ABC=[1,1,1]
        for p in f:
            h=sorted(range(4),key=lambda i:ex[i][p]); vals=sorted(ex[i][p] for i in range(4))
            low={h[0],h[1]}
            mt=0 if low in ({0,1},{2,3}) else 1 if low in ({0,2},{1,3}) else 2
            ABC[mt]*=p**(vals[2]-vals[1])
        assert ABC==[A,B,Cc],(ABC,A,B,Cc)
        x,y,z=t[0,1]*t[2,3],t[0,2]*t[1,3],t[0,3]*t[1,2]
        assert A*x-B*y+Cc*z==0
        n_=[gnorm(P) for P in Ps]
        d={key:n_[key[0]]*n_[key[1]]//g[key]**2 for key in g}
        Ds=[d[0,1]*d[2,3],d[0,2]*d[1,3],d[0,3]*d[1,2]]; Lc=lcm(*Ds)
        assert [isqrt(Lc//D) for D in Ds]==[A,B,Cc] and all(isqrt(Lc//D)**2*D==Lc for D in Ds)
        nq+=1
        eps=[data[i][2] for i in q]
        if len(set(eps))==1:
            nsame+=1
            # star angle bound |t| = sqrt(d)|sin(star angle)|, star angle = half point angle
            for key in g:
                zi,zj=data[q[key[0]]][0],data[q[key[1]]][0]
                dth=abs(math.atan2(det(zi,zj),zi[0]*zj[0]+zi[1]*zj[1]))
                maxratio=max(maxratio,abs(t[key])/(math.sqrt(d[key])*dth/2))
    print(f"{name}: n={n} k={k} C={C:.6f} reduced N={N}={dict(f)}; quartets {nq} all pass identity+matching rule+lcm formula; "
          f"same-unit {nsame}, max |t|/(sqrt(d)*angle/2) = {maxratio:.6f}")

if __name__=='__main__':
    check('8-pt record',32988780417125,8)
    check('7-pt record',235370525,7)
    check('golden n=4',537606725,8)
    check('5pt C<=2',683617125425,5)
    check('4pt C<=1',412333805,4)
    check('5pt C<=3',10414625,5)
