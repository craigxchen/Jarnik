"""[J] Sharpness of the unit constant of X.1 at g = 1: Eisenstein triples a^2 - ac + c^2 = b^2 with
q = ac squarefree are locally admissible sign-free triples with g = 1 and H = 3 q^2 (nondegenerate);
max/sqrt(q) = sqrt(max(a,c)/min(a,c)) -> 1 along m ~ 2n.  Exact check of admissibility and H."""
from math import gcd, isqrt
def sqf(n):
    p=2
    while p*p<=n:
        if n%(p*p)==0: return False
        p+=1
    return True
def primes(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0: n//=p
        p+=1
    if n>1: out.append(n)
    return out
def admissible(q,s):
    g=1
    for p in primes(q):
        z=[x%p==0 for x in s]
        if all(z): g*=p; continue
        if sum(z)!=1: return None
        o=[x*x%p for x,zz in zip(s,z) if not zz]
        if o[0]!=o[1]: return None
    return g
H=lambda a,b,c:(a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)
best=[]
for m in range(2,400):
    for n in range(1,m):
        if gcd(m,n)!=1: continue
        a,c,b=m*m-n*n,2*m*n-n*n,m*m-m*n+n*n
        if a<=0 or c<=0: continue
        g0=gcd(gcd(a,b),c); a//=g0;b//=g0;c//=g0
        q=a*c
        if q<2 or not sqf(q): continue
        g=admissible(q,(a,b,c))
        assert g is not None and H(a,b,c)==3*q*q
        best.append((max(a,c)**2/q,q,g,(a,b,c)))
best.sort()
print("count",len(best))
for r in best[:8]: print("max^2/q=%.6f q=%d g=%d sides=%s"%r)
