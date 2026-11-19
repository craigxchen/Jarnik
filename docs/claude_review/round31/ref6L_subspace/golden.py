# Record's golden 8-point family (growth/families.md Prop G): exponent check vs VI.1.
import itertools
from mpmath import mp, mpf, sqrt, log
mp.dps=80
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def ggcd(a,b):
    while b!=(0,0):
        n=b[0]**2+b[1]**2
        re=a[0]*b[0]+a[1]*b[1]; im=a[1]*b[0]-a[0]*b[1]
        q=((2*re+n)//(2*n),(2*im+n)//(2*n))
        qb=mul(q,b); a,b=b,(a[0]-qb[0],a[1]-qb[1])
    return a
def gdiv(a,b):
    n=b[0]**2+b[1]**2; c=mul(a,conj(b)); assert c[0]%n==0 and c[1]%n==0; return (c[0]//n,c[1]//n)
fib=[0,1]
for _ in range(400): fib.append(fib[-1]+fib[-2])
F=(1,2)
for n in [100,160,220]:
    J=[(fib[n+s+1],fib[n+s]) for s in range(4)]
    pts=[]
    for w in itertools.product([0,1],repeat=4):
        p=sum(w)
        if p%2==0: continue
        z=F if p==1 else conj(F)
        for s in range(4): z=mul(z, J[s] if w[s] else conj(J[s]))
        pts.append(z)
    # best 5-point subcluster by angular width; sort by argument
    import mpmath
    args=sorted(pts,key=lambda z: mpmath.atan2(z[1],z[0]))
    best=None
    for c in itertools.combinations(pts,5):
        g=c[0]
        for z in c[1:]: g=ggcd(g,z)
        prim=[gdiv(z,g) for z in c]
        H=prim[0][0]**2+prim[0][1]**2
        D=max(sqrt(mpf((u[0]-v[0])**2+(u[1]-v[1])**2)) for u in prim for v in prim)/sqrt(H)
        val=D*mpf(H)**mpf(0.25)
        if best is None or val<best[0]: best=(val,log(D)/log(H))
    print(n,"min over 5-subsets of Delta*H^(1/4) =",mpmath.nstr(best[0],12)," log Delta/log H =",mpmath.nstr(best[1],10))
