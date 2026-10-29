# Implicit equation of the round-23 hypersurface W (closure of the Pell-balanced curves C_abc) in P^4, mod p,
# by interpolation at points over F_p; then its Newton polytope, ord_e (multiplicity at e) and its cost at the
# W-vertex x_R (R = {2,3,4}: F_R = {B : |B cap R| even}), and min cost over Pi_nt. Decides WTB at the W-vertex.
import sys, itertools, random, time
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from flint import nmod_mat, nmod
import numpy as np
p=int(sys.argv[1]) if len(sys.argv)>1 else 1000003
Dmax=int(sys.argv[2]) if len(sys.argv)>2 else 14
random.seed(11)
inv=lambda v: pow(v%p,p-2,p)
def wpoint():
    while True:
        a,b,c,t=[random.randrange(1,p) for _ in range(4)]
        d=(a+b+c)%p; m=a*b*c*d%p
        sq=[a*a%p,b*b%p,c*c%p,d*d%p]; S=sum(sq)%p
        e3=sum(sq[i]*sq[j]*sq[k] for i,j,k in itertools.combinations(range(4),3))%p
        B=(e3-S*m)*inv(4)%p
        za=(a*a*b+a*a*c+a*b*b+a*c*c-b*b*c-b*c*c)*inv(2)%p
        assert (za*za-(m*a*a+B))%p==0
        den=(t*t-m)%p
        if den==0: continue
        x=(a*(t*t+m)-2*za*t)*inv(den)%p
        z=(za+t*(x-a))%p
        assert (z*z-(m*x*x+B))%p==0
        F0=(x-a)*(x-b)*(x-c)*(x+d)%p
        F1=(pow(x,4,p)-S*inv(2)*x*x+m+2*x*z)%p
        F2=(x*x-b*b)*(x*x-c*c)%p; F3=(x*x-a*a)*(x*x-c*c)%p; F4=(x*x-a*a)*(x*x-b*b)%p
        F=[F0,F1,F2,F3,F4]
        # check norm identity F1*F1bar = n
        F1b=(pow(x,4,p)-S*inv(2)*x*x+m-2*x*z)%p
        n=(x*x-a*a)*(x*x-b*b)*(x*x-c*c)*(x*x-d*d)%p
        assert F1*F1b%p==n
        if all(F): return F
def monos(D): return [e for e in itertools.product(range(D+1),repeat=5) if sum(e)==D]
for D in range(2,Dmax+1):
    ms=monos(D); N=len(ms)
    pts=[wpoint() for _ in range(N+40)]
    t0=time.time()
    rows=[[ (pow(P[0],e[0],p)*pow(P[1],e[1],p)%p*pow(P[2],e[2],p)%p*pow(P[3],e[3],p)%p*pow(P[4],e[4],p))%p for e in ms] for P in pts]
    M=nmod_mat(rows,p)
    r=M.rank()
    print("D=",D,"monomials",N,"rank",r,"kernel dim",N-r,f"({time.time()-t0:.1f}s)",flush=True)
    if r<N:
        # kernel vector: nullspace of M (rows=points) -> columns combos
        K,nul=M.nullspace()
        vec=[int(K[i,0]) for i in range(N)]
        supp=[ms[i] for i in range(N) if vec[i]%p!=0]
        print("support size",len(supp))
        np.save('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/w_support_D%d_p%d.npy'%(D,p),np.array(supp))
        np.save('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/w_coeffs_D%d_p%d.npy'%(D,p),np.array([vec[i] for i in range(N) if vec[i]%p!=0]))
        break
