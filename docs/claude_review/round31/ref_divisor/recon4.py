# Reconstruction test on a genuine 4-point cluster (construction copied from Astra's near_square_check, analysis independent).
import sys
from functools import reduce
from math import isqrt, gcd
from itertools import combinations
from fractions import Fraction
from mpmath import mp, atan2, mpf
mp.dps=600
def mul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0],-a[1])
def nrm(a): return a[0]**2+a[1]**2
def ggcd(z,w):
    while w!=(0,0):
        n=nrm(w); x,y=mul(z,cj(w)); q=((2*x+n)//(2*n),(2*y+n)//(2*n)); a,b=mul(q,w); z,w=w,(z[0]-a,z[1]-b)
    return z
def div(z,w):
    n=nrm(w); x,y=mul(z,cj(w)); assert x%n==0 and y%n==0; return (x//n,y//n)
n=10**30+1; x,y,z,u=-3,-2,6,-6
blocks=[]
for a,b,c in [(x,y,z),(x,z,y),(y,z,x)]:
    blocks.append({(0,1),(u,a-1),(u,b-1),(u*(1+c)//c,c)})
roots=sorted(set.union(*blocks)); pts=[]
for ch in [set(),*blocks]:
    f=[(n-4*re,(-4 if (re,im) in ch else 4)*im) for re,im in roots]
    pts.append(reduce(mul,f,(1,0)))
g=reduce(ggcd,pts); Z=[div(p,g) for p in pts]; N=nrm(Z[0]); assert all(nrm(w)==N for w in Z)
m=len(Z)
# edge data from divisors only
A={}
for i,j in combinations(range(m),2):
    d=N//nrm(ggcd(Z[i],Z[j])); s=isqrt(d); t2=d-s*s; t=isqrt(t2)
    assert t>0 and t*t==t2 and gcd(s,t)==1 and 256*t**4*N<=d*d, (i,j)
    A[(i,j)]=(s,t)
print("all 6 edges pass primitive floor test; labels distinct:",len({nrm(v) for v in A.values()})==6)
psi={e:atan2(v[1],v[0]) for e,v in A.items()}
def P(i,j): return psi[(min(i,j),max(i,j))]
# triangle theorem check exact
for i,j,k in combinations(range(m),3):
    a,b,c=A[(i,j)],A[(j,k)],A[(i,k)]
    ok=any(mul(mul(p,q),cj(r))[1]==0 and mul(mul(p,q),cj(r))[0]>0 for p,q,r in ((a,b,c),(a,c,b),(b,c,a)))
    print("triangle",(i,j,k),"exact degenerate:",ok)
D=max(psi.values()); diam=[e for e,v in psi.items() if v==D]; print("diameter edges",diam)
Aend,Bend=diam[0]
xs={i:(0 if i==Aend else P(Aend,i)) for i in range(m)}
# line embedding check
okline=all(abs(P(i,j)-abs(xs[i]-xs[j]))<mpf(10)**-500 for i,j in combinations(range(m),2))
print("line embedding consistent:",okline)
# w_i = e^{2 i x_i} in Q(i): w_i = A/conj(A) for edge (Aend,i) with orientation from sign
recon=[]
for i in range(m):
    if i==Aend: recon.append(Z[Aend]); continue
    e=(min(Aend,i),max(Aend,i)); a=A[e]
    # candidate: z_i = z_A * (a/conj a)^{+-1}
    found=None
    for num,den in ((a,cj(a)),(cj(a),a)):
        X=mul(Z[Aend],num); nn=nrm(den); q=mul(X,cj(den))
        if q[0]%nn==0 and q[1]%nn==0:
            c=(q[0]//nn,q[1]//nn)
            for uu in ((1,0),(0,1),(-1,0),(0,-1)):
                if mul(c,uu)==Z[i]: found=(num==a,uu)
    recon.append(found)
print("reconstructed points equal actual up to unit:",recon[1:] if Aend==0 else recon)
