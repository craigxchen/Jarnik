"""Sanity checks (evidence; the proofs do not use them):
(a) Prop 3.3(1)-(2) and Lemma 3.2 on genuine Gaussian prime angles: max arc counts vs the lattice bound, and the
    optimal large-sieve constant (largest eigenvalue of the Dirichlet-kernel Gram matrix) vs 8Km and the converse.
(b) Prop 4.2 demo: a random 0/1 profile (M rows, s columns, window primes); random points of the realisation torus
    (phi0 + H°) satisfy the realisation equations and are spread."""
import math, numpy as np
from sympy import primerange, Matrix
from glib import gprime
rng=np.random.default_rng(7)
Q=math.pi/2
def arcmax(ang, eta):
    a=np.sort(np.mod(ang,Q)); n=len(a); b=np.concatenate([a,a+Q]); j=0; best=0
    for i in range(n):
        while j<i+n and b[j]<a[i]+eta: j+=1
        best=max(best,j-i)
    return best
def lsconst(ang,K):
    d=ang[:,None]-ang[None,:]
    G=np.ones_like(d)
    for k in range(1,K+1): G+=2*np.cos(4*k*d)
    return np.linalg.eigvalsh(G)[-1]
print("(a) genuine primes")
for Y in [20000,100000]:
    ps=[p for p in primerange(5,Y+1) if p%4==1]
    th=np.array([math.atan2(*gprime(p)[::-1]) for p in ps])
    sub=th if len(th)<=2500 else rng.choice(th,2500,replace=False)
    for K in [10,int(Y**0.5),300]:
        eta=math.pi/(4*K)
        m=arcmax(th,eta); bound=eta*Y/2+2.53*Y**0.5+1.58
        lam=lsconst(sub,K); msub=arcmax(sub,eta); t=arcmax(sub,eta/2)
        print(" Y=%d n=%d K=%d: maxarc=%d <= %.1f ; subset n=%d: Lambda=%.1f <= 8Km=%.1f, >= (2K+1)t/sqrt2=%.1f; genuine bound piY+20.24K sqrtY+12.64K=%.1f"%(
            Y,len(th),K,m,bound,len(sub),lam,8*K*msub,(2*K+1)*t/math.sqrt(2),math.pi*Y+20.24*K*Y**0.5+12.64*K))
print("(b) realisation torus demo")
M=40; ps=[p for p in primerange(int(M**1.3),int(M**1.9)) if p%4==1][:60]; s=len(ps)
A=rng.integers(0,2,size=(M,s))
Bn=2.0*(A[1:]-A[0])
_,sv,Vt=np.linalg.svd(Bn); rk=int((sv>1e-9).sum()); U=Vt[rk:]   # real orthonormal basis of ker B (a rational subspace)
# Haar on H° approximated by phi0 + U^T t with t uniform in a large box (equidistribution of the rational subtorus)
class _R: pass
B=_R(); B.rank=lambda: rk
print(" M=%d s=%d rank B=%d dim H=%d primes in [%d,%d]"%(M,s,B.rank(),len(U),ps[0],ps[-1]))
e=np.ones(s)  # e_j = 1 (squarefree)
phi0=rng.random(s)*1e-9
for trial in range(3):
    t=rng.random(len(U))*1000.0
    phi=phi0+t@U
    ang=np.mod((2*A-e)@phi,Q); spread=np.ptp(np.mod(ang-ang[0]+Q/4,Q))
    out=[]
    for L in [2,4,6]:
        eta=Q*2.0**-L; out.append((L,arcmax(phi,eta),round((16/math.pi)*s*eta+4*math.sqrt(L*s),1)))
    print(" trial %d: realisation spread of the M angles mod pi/2 = %.2e ; (L, max arc count, Prop 4.2 bound):"%(trial,spread),out)
