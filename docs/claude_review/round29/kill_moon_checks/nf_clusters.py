"""Independent recomputation (kill_moon): record clusters, reduced chords, odd-prime collision counts,
Erdos-Turan usefulness test n < M^2/88.7, power-sum product formula (with the factor 4)."""
import math, itertools, cmath
from sympy import factorint, isprime
from fractions import Fraction
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
def nm(a): return a[0]**2+a[1]**2
def gg(a,b):
    while b!=(0,0):
        n=nm(b); t=mul(a,cj(b))
        q=(round(Fraction(t[0],n)),round(Fraction(t[1],n))); r=mul(b,q)
        a,b=b,(a[0]-r[0],a[1]-r[1])
    return a
def gdiv(a,b):
    n=nm(b); t=mul(a,cj(b)); assert t[0]%n==0 and t[1]%n==0; return (t[0]//n,t[1]//n)
def cornacchia(p):
    # p=1 mod 4: find x with x^2=-1 mod p then Euclid
    for g in range(2,p):
        x=pow(g,(p-1)//4,p)
        if x*x%p==p-1: break
    a,b=p,x
    while b*b>p: a,b=b,a%b
    c=math.isqrt(p-b*b); assert b*b+c*c==p; return (b,c)
def allpoints(N):
    pts=[(1,0)]
    for p,e in factorint(N).items():
        if p==2: 
            for _ in range(e): pts=[mul(z,(1,1)) for z in pts]
        elif p%4==3: assert e%2==0; pts=[(z[0]*p**(e//2),z[1]*p**(e//2)) for z in pts]
        else:
            pi=cornacchia(p); new=[]
            for z in pts:
                for a in range(e+1):
                    w=z
                    for _ in range(a): w=mul(w,pi)
                    for _ in range(e-a): w=mul(w,cj(pi))
                    new.append(w)
            pts=new
    out=set()
    for z in pts:
        w=z
        for _ in range(4): out.add(w); w=mul(w,(0,1))
    for z in out: assert nm(z)==N
    return sorted(out,key=lambda z: math.atan2(z[1],z[0]))
for N,k in [(537606725,8),(32988780417125,8),(348984761607530165,8),(235370525,7),(683617125425,5)]:
    P=allpoints(N); n=len(P); R=math.sqrt(N)
    ang=[math.atan2(z[1],z[0]) for z in P]
    best=min(((ang[(s+k-1)%n]-ang[s])%(2*math.pi),s) for s in range(n))
    cl=[P[(best[1]+j)%n] for j in range(k)]
    C=best[0]*R/math.sqrt(R)
    # chord-based C (max chord / sqrt R)
    Cch=max(math.sqrt(nm((x[0]-y[0],x[1]-y[1]))) for x in cl for y in cl)/math.sqrt(R)
    g=cl[0]
    for z in cl[1:]: g=gg(g,z)
    cl2=[gdiv(z,g) for z in cl]
    coll={}
    big={}
    for x,y in itertools.combinations(cl2,2):
        d=(x[0]-y[0],x[1]-y[1]); c=gdiv(d,gg(x,y))
        for l,e in factorint(nm(c)).items():
            if l==2: continue
            if l<=4*k: coll[l]=coll.get(l,0)+1
            else: big[l]=big.get(l,0)+1
    print(f"N={N} n=r(N)={n} k={k} C_arc={C:.6f} C_chord={Cch:.6f} gcd norm={nm(g)}")
    print("   odd l<=4k: pairs",dict(sorted(coll.items())),"; max pairs at a prime l>4k:",max(big.values()) if big else 0)
    print(f"   ET useful only if n < M^2/88.7: {n} < {k*k/88.7:.2f}? {n<k*k/88.7}")
    # power sum check m=1..3 vs product formula with factor 4 (on gcd-normalised N? use full N with 2 and inert handled)
    th=[math.atan2(z[1],z[0]) for z in P]
    f=factorint(N); 
    for m in (1,2,3):
        s=sum(cmath.exp(4j*m*t) for t in th)
        prod=4
        for p,e in f.items():
            if p%4==1:
                pi=cornacchia(p); phi=math.atan2(pi[1],pi[0]); x=4*m*phi
                prod*= sum(cmath.exp(1j*x*(2*a-e)) for a in range(e+1))
        # 2-part rotates all points by (1+i)^e2 -> exp(4im*e2*pi/4)=(-1)^{m e2}
        e2=f.get(2,0); prod*=(-1)**(m*e2)
        assert abs(s-prod)<1e-6*n, (s,prod)
    print("   power sums s_m/R^{4m} = 4*prod_p sin((e+1)x)/sin(x) verified for m=1,2,3 (factor 4 required)")
