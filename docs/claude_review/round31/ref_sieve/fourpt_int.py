# Explicit integer instances of the section-6 family at x=-3,m=3 (d=4), check primary, gcd, exponent, t_ab
from mpmath import mp, mpf, mpc, arg, sqrt, log
mp.dps=80
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivround(a,b):
    n=gnorm(b); num=gmul(a,(b[0],-b[1]))
    def rd(p): return (2*p+n)//(2*n)
    return (rd(num[0]),rd(num[1]))
def gsub(a,b): return (a[0]-b[0],a[1]-b[1])
def ggcd(a,b):
    while b!=(0,0):
        q=gdivround(a,b); a,b=b,gsub(a,gmul(q,b))
    return a
def gexact(a,b):
    n=gnorm(b); num=gmul(a,(b[0],-b[1])); assert num[0]%n==0 and num[1]%n==0; return (num[0]//n,num[1]//n)
def primary_assoc(a):
    for u in [(1,0),(0,1),(-1,0),(0,-1)]:
        b=gmul(a,u)
        # primary: b = 1 mod (2+2i)  <=> b-1 divisible by 2+2i
        c=gsub(b,(1,0)); n=8; num=gmul(c,(2,-2))
        if num[0]%8==0 and num[1]%8==0: return b
    raise Exception('even')
P1=[(0,1),(-6,-4),(-6,-3),(-7,6)]      # (a,b,c)=(x,y,z): roots i, u+i(x-1), u+i(y-1), u(1+z)/z + i z
P2=[(0,1),(-6,-4),(-6,5),(-3,-2)]      # (x,z,y)
P3=[(0,1),(-6,-3),(-6,5),(-4,-3)]      # (y,z,x)
roots=[(0,1),(-6,-4),(-6,-3),(-6,5),(-7,6),(-3,-2),(-4,-3)]
for P in (P1,P2,P3): assert set(P)<=set(roots)
d=4
for n in [10**3+1, 10**5+1, 10**7+1, 10**9+1, 10**12+1]:
    H={r:(n-d*r[0],-d*r[1]) for r in roots}
    Hb={r:(H[r][0],-H[r][1]) for r in roots}
    pts=[]
    z0=(1,0)
    for r in roots: z0=gmul(z0,Hb[r])
    pts.append(z0)
    for P in (P1,P2,P3):
        zz=(1,0)
        for r in roots: zz=gmul(zz, H[r] if r in P else Hb[r])
        pts.append(zz)
    g=pts[0]
    for q in pts[1:]: g=ggcd(g,q)
    g=primary_assoc(g)
    npts=[gexact(q,g) for q in pts]
    assert all(primary_assoc(q)==q for q in npts)
    N=gnorm(npts[0]); assert all(gnorm(q)==N for q in npts)
    R=sqrt(mpf(N))
    zc=[mpc(q[0],q[1]) for q in npts]
    diam=max(abs(a-b) for a in zc for b in zc)
    ts=[]
    for i in range(4):
        for j in range(i+1,4):
            gij=ggcd(npts[i],npts[j]); dij=N//gnorm(gij)
            th=arg(zc[i]/zc[j]); t=sqrt(mpf(dij))*mp.sin(th/2)
            ts.append(t)
    print("  diam/n^3=",mp.nstr(diam/mpf(n)**3,8)," R/n^7=",mp.nstr(R/mpf(n)**7,8)); print(f"n={n}: N(g)={gnorm(g)}, log R/log n={float(log(R)/log(n)):.4f}, log diam/log R={float(log(diam)/log(R)):.4f}, t_ab=",[mp.nstr(t,8) for t in ts])
