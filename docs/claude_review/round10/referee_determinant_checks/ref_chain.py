# Referee: exact end-to-end check of the Liouville/forced-divisibility chain of Theorem 4.1 on actual
# lattice clusters (exact Gaussian integers; mpmath only for the archimedean side).
# For integer c in V(A,m) (exact rational nullspace, scaled to integers):
#  (L) if F_c(z) != 0:  Norm F_c(z) * prod_p p^{h_{supp c}(a_p)} >= N^D        (exact integer test)
#  (U) |F_c(z)| <= N^{D/2} sum_k |c_k| |<k,delta>|^m / m!  <= N^{D/2} |c|_1 (D Theta/2)^m / m!
#  (P) F_c(z) == N^{D/2} sum_k c_k e^{i<k,delta>}   (phase convention; delta measured from arc midpoint)
import itertools, random, math, sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 60
random.seed(2026)

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gpow(a,e):
    r=(1,0); b=a
    while e:
        if e&1: r=gmul(r,b)
        b=gmul(b,b); e>>=1
    return r
def gdivisible(a,b):
    n=b[0]**2+b[1]**2; t=gmul(a,gconj(b)); return t[0]%n==0 and t[1]%n==0
def gdivexact(a,b):
    n=b[0]**2+b[1]**2; t=gmul(a,gconj(b)); return (t[0]//n,t[1]//n)
def vpi(a,pi):
    v=0
    while gdivisible(a,pi): a=gdivexact(a,pi); v+=1
    return v
def gprime(p):
    for x in range(1,p):
        y=math.isqrt(p-x*x)
        if y*y==p-x*x: return (x,y)

def lattice_points(N):
    pts=[]
    for x in range(-math.isqrt(N), math.isqrt(N)+1):
        y2=N-x*x; y=math.isqrt(y2)
        if y*y==y2:
            pts.append((x,y));
            if y>0: pts.append((x,-y))
    pts=sorted(set(pts), key=lambda z: math.atan2(z[1],z[0]))
    return pts

out=[]
def pr(s): print(s); out.append(s)

cases=[ ({5:1,13:1,17:1,29:1},4,2,3), ({5:2,13:1,17:1,29:1},4,2,3), ({5:1,13:1,17:1,29:1,37:1},5,1,3),
        ({5:2,13:2,17:1},4,2,4), ({5:1,13:1,17:1,29:1,37:1,41:1},5,2,3) ]
tot=dict(L=0,U=0,P=0,zero=0,nonzero=0,fail=0)
for fac,n,T,m in cases:
    N=1
    for p,e in fac.items(): N*=p**e
    pts=lattice_points(N)
    pis={p:gprime(p) for p in fac}
    for trial in range(12):
        i0=random.randrange(len(pts)-n); Z=pts[i0:i0+n]
        th=[mp.atan2(z[1],z[0]) for z in Z]
        th0=(max(th)+min(th))/2; dl=[t-th0 for t in th]; Theta=max(th)-min(th)
        box=[tuple(list(k)+[-sum(k)]) for k in itertools.product(range(-T,T+1),repeat=n-1)]
        A=random.sample(box, min(len(box), 45))
        mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=n-1) if sum(a)==t]
        Mx=sp.Matrix([[sp.prod([k[j]**a[j] for j in range(n-1)]) for k in A] for a in mons])
        ns=Mx.nullspace()
        D=max(sum(abs(t) for t in k) for k in A); D+=D%2
        for vec in ns[:4]:
            den=sp.ilcm(*[sp.fraction(x)[1] for x in vec]); c=[int(x*den) for x in vec]
            if all(x==0 for x in c): continue
            F=(0,0)
            for ck,k in zip(c,A):
                if ck==0: continue
                g=gpow((N,0),(D-sum(abs(t) for t in k))//2)
                for x,kx in enumerate(k):
                    g=gmul(g, gpow(Z[x],kx) if kx>0 else gpow(gconj(Z[x]),-kx))
                F=(F[0]+ck*g[0],F[1]+ck*g[1])
            supp=[k for ck,k in zip(c,A) if ck!=0]
            # (P) phase identity
            S=mp.mpc(0)
            for ck,k in zip(c,A):
                if ck: S+=ck*mp.e**(1j*mp.fsum(kx*d for kx,d in zip(k,dl)))
            Fm=mp.mpc(F[0],F[1]); scale=mp.mpf(N)**(mp.mpf(D)/2)
            okP = abs(Fm-scale*S) <= mp.mpf(10)**(-30)*scale*(1+sum(abs(x) for x in c))
            tot['P']+=okP
            # (U) Taylor
            ub1=scale*mp.fsum(abs(ck)*abs(mp.fsum(kx*d for kx,d in zip(k,dl)))**m for ck,k in zip(c,A))/mp.factorial(m)
            ub2=scale*sum(abs(x) for x in c)*(D*Theta/2)**m/mp.factorial(m)
            okU = abs(Fm) <= ub1*(1+mp.mpf(10)**-40) and ub1 <= ub2*(1+mp.mpf(10)**-40)
            tot['U']+=okU
            # (L) forced divisibility
            if F==(0,0):
                tot['zero']+=1; okL=True
            else:
                tot['nonzero']+=1
                Nm=F[0]**2+F[1]**2; prodp=1
                for p,e in fac.items():
                    a=[vpi(z,pis[p]) for z in Z]
                    vals=[sum(kx*ax for kx,ax in zip(k,a)) for k in supp]
                    prodp*=p**(max(vals)-min(vals))
                okL = Nm*prodp >= N**D
            tot['L']+=okL
            if not (okP and okU and okL):
                tot['fail']+=1; pr("FAIL %s %s %s %s"%(fac,Z,okP,(okU,okL)))
    pr("case N=%s n=%d T=%d m=%d done"%(fac,n,T,m))
pr("totals: %s"%tot)
open(sys.argv[1],"w").write("\n".join(out)+"\n")
