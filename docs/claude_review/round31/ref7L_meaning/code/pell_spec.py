# Independent check (own arithmetic, no sympy): polynomial Pell five-mark family of Part VII.5,
# (a) D,h over Q(i)[t]; (b) integer specialisation t=T: primitive Gaussian 5-tuples, exponent log Delta/log H;
# (c) Jacobian rank of seed -> anchored arguments map at fixed T (density of the specialised family).
from fractions import Fraction as F
import math, sys, itertools
# Gaussian rationals as tuples (re,im) of Fractions
def gadd(a,b): return (a[0]+b[0],a[1]+b[1])
def gsub(a,b): return (a[0]-b[0],a[1]-b[1])
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def ginv(a):
    n=a[0]*a[0]+a[1]*a[1]; return (a[0]/n,-a[1]/n)
Z=(F(0),F(0)); ONE=(F(1),F(0)); I=(F(0),F(1))
def norm(p):
    while p and p[-1]==Z: p=p[:-1]
    return p
def padd(p,q):
    n=max(len(p),len(q)); return norm([gadd(p[i] if i<len(p) else Z, q[i] if i<len(q) else Z) for i in range(n)])
def psub(p,q):
    n=max(len(p),len(q)); return norm([gsub(p[i] if i<len(p) else Z, q[i] if i<len(q) else Z) for i in range(n)])
def pmul(p,q):
    if not p or not q: return []
    r=[Z]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        if a==Z: continue
        for j,b in enumerate(q): r[i+j]=gadd(r[i+j],gmul(a,b))
    return norm(r)
def pconj(p): return [gconj(a) for a in p]
def pdivmod(a,b):
    a=list(a); q=[Z]*max(1,len(a)-len(b)+1); ib=ginv(b[-1])
    while len(a)>=len(b) and a:
        c=gmul(a[-1],ib); k=len(a)-len(b); q[k]=c
        for i,bb in enumerate(b): a[i+k]=gsub(a[i+k],gmul(c,bb))
        a=norm(a[:-1]) if a[-1]==Z else norm(a)
    return norm(q),a
def pgcd(a,b):
    while b: a,b=b,pdivmod(a,b)[1]
    return [gmul(c,ginv(a[-1])) for c in a]
def deg(p): return len(p)-1
def C(x): return (F(x),F(0))
t=[Z,ONE]
def family(n,seeds,edges):
    A=[[padd(pmul(t,t),[ONE]),t],[t,[ONE]]]
    G=[]
    for s in seeds:
        x,y=[ONE],[C(s)]
        for _ in range(n):
            x,y=padd(pmul(A[0][0],x),pmul(A[0][1],y)),padd(pmul(A[1][0],x),pmul(A[1][1],y))
        G.append(padd(x,pmul([I],y)))
    beta=[(F(0),F(2)),ONE]
    z0=beta
    for g in G: z0=pmul(z0,pconj(g))
    marks=[z0]
    for (j,k) in edges:
        z=pconj(beta)
        for l,g in enumerate(G): z=pmul(z, g if l in (j,k) else pconj(g))
        marks.append(z)
    return marks
def check_poly(n,seeds,edges):
    m=family(n,seeds,edges)
    N=pmul(m[0],pconj(m[0]))
    assert all(pmul(z,pconj(z))==N for z in m)
    g=m[0]
    for z in m[1:]: g=pgcd(g,z)
    D=deg(m[0]); assert all(deg(z)==D for z in m)
    dd=max(deg(psub(a,b)) for a,b in itertools.combinations(m,2))
    return D,deg(g),D-dd
# Gaussian integers
def ev(p,T):  # exact evaluation, coefficients Gaussian rationals with denominator 1 here
    re=F(0); im=F(0); pw=F(1)
    for c in p: re+=c[0]*pw; im+=c[1]*pw; pw*=T
    assert re.denominator==1 and im.denominator==1
    return (int(re),int(im))
def evf(p,T):
    re=F(0); im=F(0); pw=F(1)
    for c in p: re+=c[0]*pw; im+=c[1]*pw; pw*=T
    return complex(float(re),float(im))
def zmul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def zdivround(a,b):
    n=b[0]*b[0]+b[1]*b[1]; num=zmul(a,(b[0],-b[1]))
    def rd(x): return (2*x+n)//(2*n)
    return (rd(num[0]),rd(num[1]))
def zmod(a,b):
    q=zdivround(a,b); qb=zmul(q,b); return (a[0]-qb[0],a[1]-qb[1])
def zgcd(a,b):
    while b!=(0,0): a,b=b,zmod(a,b)
    return a
def zexact(a,b):
    q=zdivround(a,b); assert zmul(q,b)==a; return q
def specialise(n,seeds,edges,T):
    m=family(n,seeds,edges); pts=[ev(z,T) for z in m]
    g=pts[0]
    for p in pts[1:]: g=zgcd(g,p)
    pts=[zexact(p,g) for p in pts]
    H=pts[0][0]**2+pts[0][1]**2
    assert all(p[0]**2+p[1]**2==H for p in pts) and len(set(pts))==5
    lH=math.log(H)
    md=max(math.log(math.hypot(a[0]-b[0],a[1]-b[1])) for a,b in itertools.combinations(pts,2))
    logDelta=md-lH/2
    return g[0]**2+g[1]**2, lH/math.log(10), logDelta/lH
def angles(n,seeds,edges,T):
    import cmath
    m=family(n,seeds,edges)
    # float evaluation via exact ints then angle of ratio to anchor
    pts=[evf(z,T) for z in m]
    a0=pts[0]
    return [cmath.phase(p/a0) for p in pts[1:]]
def jac_rank(n,seeds,edges,T,h=F(1,10**6)):
    # finite differences of the anchored arguments with respect to rational seeds (exact pts, float angles)
    import numpy
    base=angles(n,seeds,edges,T); J=[]
    for i in range(len(seeds)):
        s2=list(seeds); s2[i]=s2[i]+h
        J.append([(x-y)/float(h) for x,y in zip(angles(n,s2,edges,T),base)])
    sv=numpy.linalg.svd(numpy.array(J),compute_uv=False)
    return sv
if __name__=="__main__":
    E5=[(0,1),(0,2),(1,2),(0,3)]
    seeds=[1,3,5,8]
    for n in [1,2,3]:
        print("poly n=%d (D,deg gcd,h)="%n, check_poly(n,seeds,E5), "expected",(8*n+1,0,4*n))
    for n in [1,2,3]:
        for T in [10,100,1000,10**4]:
            Ng,lgH,ex=specialise(n,seeds,E5,T)
            print("n=%d T=%d N(gcd)=%d log10H=%.1f logDelta/logH=%.4f  target -h/(2D)=%.4f"%(n,T,Ng,lgH,ex,-4*n/(2*(8*n+1))))
