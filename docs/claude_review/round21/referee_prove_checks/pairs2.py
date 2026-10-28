# Independent recomputation of prove.md §7.2: record clusters, profiles, rank B, Smith invariants, D_j.
import math, itertools
from sympy import Matrix, factorint
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
from mpmath import mp, mpf, atan2, pi as mpi
mp.dps = 50

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivexact(a,b):
    # a/b if exact else None
    n=gnorm(b); num=gmul(a,gconj(b))
    if num[0]%n or num[1]%n: return None
    return (num[0]//n, num[1]//n)
def gaussprime(p):
    for x in range(1,int(math.isqrt(p))+1):
        y2=p-x*x; y=math.isqrt(y2)
        if y*y==y2 and y>0: return (x,y)
def points(n):
    pts=[]
    for x in range(0, math.isqrt(n)+1):
        y2=n-x*x; y=math.isqrt(y2)
        if y*y==y2: pts.append((x,y))
    # all four quadrants
    allp=set()
    for (x,y) in pts:
        for s in [(x,y),(-x,y),(x,-y),(-x,-y)]: allp.add(s)
    return sorted(allp, key=lambda z: math.atan2(z[1],z[0]))
def best_window(n,k):
    P=points(n); ang=[mp.atan2(z[1],z[0]) for z in P]; m=len(P)
    R=mp.sqrt(n); best=None
    for i in range(m):
        j=(i+k-1)%m
        d=ang[j]-ang[i]
        if d<0: d+=2*mpi
        if best is None or d<best[0]: best=(d,i)
    d,i=best
    win=[P[(i+t)%m] for t in range(k)]
    C_arc=R*d/mp.sqrt(R)
    chord=mp.sqrt(gnorm((win[-1][0]-win[0][0],win[-1][1]-win[0][1])))
    return win, C_arc, chord/mp.sqrt(R)
def val(u,pi_):
    a=0
    while True:
        q=gdivexact(u,pi_)
        if q is None: return a
        u=q; a+=1
def ggcd(a,b):
    while b!=(0,0):
        n=gnorm(b); num=gmul(a,gconj(b))
        q=(round(num[0]/n), round(num[1]/n))
        r=(a[0]-(b[0]*q[0]-b[1]*q[1]), a[1]-(b[0]*q[1]+b[1]*q[0]))
        a,b=b,r
    return a
def hnf_last(Bm, j):
    # order of e_j in Z^s / rowlattice(B): permute col j to last, row-echelon over Z
    s=len(Bm[0]); cols=[c for c in range(s) if c!=j]+[j]
    A=[[row[c] for c in cols] for row in Bm]
    r=0
    for c in range(s):
        # gcd-reduce column c among rows r..
        while True:
            nz=[i for i in range(r,len(A)) if A[i][c]!=0]
            if len(nz)<=1: break
            i0=min(nz,key=lambda i:abs(A[i][c]))
            for i in nz:
                if i!=i0:
                    q=A[i][c]//A[i0][c]
                    A[i]=[A[i][t]-q*A[i0][t] for t in range(s)]
        nz=[i for i in range(r,len(A)) if A[i][c]!=0]
        if not nz: continue
        A[r],A[nz[0]]=A[nz[0]],A[r]; r+=1
    # rows with pivot in last column
    for row in A:
        if all(x==0 for x in row[:-1]) and row[-1]!=0: return abs(row[-1])
    return math.inf
for (name,n,k) in [("golden8",537606725,8),("5pt",2163838625,5)]:
    win,Carc,Cch=best_window(n,k)
    G=win[0]
    for z in win[1:]: G=ggcd(G,z)
    Np=n//gnorm(G)
    u=[gdivexact(z,G) for z in win]
    f=factorint(Np); primes=sorted(f)
    pis=[gaussprime(p) for p in primes]
    A=[[val(z,pi_) for pi_ in pis] for z in u]
    Bm=[[2*(A[x][j]-A[0][j]) for j in range(len(primes))] for x in range(1,k)]
    rk=Matrix(Bm).rank()
    snf=smith_normal_form(Matrix(Bm).T, domain=ZZ)
    inv=[abs(snf[i,i]) for i in range(min(snf.shape))]
    D=[hnf_last(Bm,j) for j in range(len(primes))]
    sp=sum(p**-0.5 for p in primes)
    print(name,"n=",n,"k=",k,"C_arc=%.4f C_chord=%.4f"%(float(Carc),float(Cch)),"G=",G,"N'=",f,"s=",len(primes),"P+=",max(primes),
          "P+/k^2=%.3f"%(max(primes)/k**2),"sum p^-1/2=%.4f"%sp,"rankB=",rk,"Smith(B^T)=",inv,"D_j=",D)
    s=len(primes); dep=[]
    for i,j in itertools.combinations(range(s),2):
        keep=[c for c in range(s) if c not in (i,j)]
        r2=Matrix([[row[c] for c in keep] for row in Bm]).rank() if keep else 0
        if rk-r2>0: dep.append((primes[i],primes[j]))
    print(name,"dependent pairs:",len(dep),"of",s*(s-1)//2,dep)
