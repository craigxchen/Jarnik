# doubled k=4 V-null (4,4): columns 1,2,3 (+5 normalised); each B subset {1,2,3} carries two places rho_B^(1), rho_B^(2)
import itertools, numpy as np, sys
from scipy.optimize import least_squares
B=[frozenset(c) for k in range(1,4) for c in itertools.combinations(range(1,4),k)]  # 7
# unknown places: index 2*i, 2*i+1 for B[i]
MEM=[np.array([2*i+t for i,b in enumerate(B) if j in b for t in (0,1)]) for j in range(1,4)]
K=np.array([1,2,3,5,6,7,8])
def polyjac(roots):
    c=np.poly(roots); n=len(roots); J=np.zeros((n+1,n),dtype=complex)
    for i in range(n): J[1:,i]=-np.poly(np.delete(roots,i))
    return c,J
class P:
    def __init__(s,eps,rng,nsl):
        s.eps=np.array(eps,float); s.L=rng.normal(size=(nsl,28)); s.c0=rng.normal(size=nsl)
    def rho(s,x): return x[:14]+1j*s.eps*np.exp(x[14:])
    def F(s,x):
        r=s.rho(x); out=[np.poly(r[MEM[j]])[K].imag for j in range(3)]
        return np.concatenate(out+[[np.log(abs(r[13]))],0.1*(s.L@x-s.c0)])
    def J(s,x):
        r=s.rho(x); da=np.ones(14,complex); db=1j*s.eps*np.exp(x[14:]); rows=[]
        for j in range(3):
            c,Jc=polyjac(r[MEM[j]]); blk=np.zeros((7,28)); sub=Jc[K]
            blk[:,MEM[j]]=(sub*da[MEM[j]]).imag; blk[:,14+MEM[j]]=(sub*db[MEM[j]]).imag; rows.append(blk)
        g=np.zeros((1,28)); g[0,13]=(np.conj(r[13])*da[13]).real/abs(r[13])**2; g[0,27]=(np.conj(r[13])*db[13]).real/abs(r[13])**2
        return np.vstack(rows+[g,0.1*s.L])
def eqJ(r):  # jacobian of the 21 equations wrt (Re rho, Im rho) for rank computation
    rows=[]
    for j in range(3):
        c,Jc=polyjac(r[MEM[j]]); blk=np.zeros((7,28)); sub=Jc[K]
        blk[:,MEM[j]]=sub.imag; blk[:,14+MEM[j]]=(1j*sub).imag; rows.append(blk)
    return np.vstack(rows)
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
    # sign patterns: each R_j 4 up / 4 down
    pats=[e for e in itertools.product([1,-1],repeat=14) if all(sum(e[m] for m in MEM[j])==0 for j in range(3))]
    nd=0
    for it in range(N):
        eps=pats[rng.integers(len(pats))]
        pr=P(eps,rng,28-21-1)
        x0=np.concatenate([rng.normal(size=14)*rng.choice([0.5,1,2]),rng.normal(size=14)*0.7])
        try: sol=least_squares(pr.F,x0,jac=pr.J,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=1500)
        except Exception: continue
        r=pr.rho(sol.x); res=np.linalg.norm(pr.F(sol.x)[:21])
        allv=np.concatenate([r,np.conj(r)]); iu=np.triu_indices(28,1); sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu])
        if res>1e-10 or np.min(np.abs(r.imag))<1e-3 or sep<1e-3 or np.max(np.abs(r))>1e4: continue
        nd+=1
        sv=np.linalg.svd(eqJ(r),compute_uv=False); rank=int(np.sum(sv>1e-8*sv[0]))
        pull=max(min(abs(r[2*i]+r[2*i+1]), abs(r[2*i]+np.conj(r[2*i+1])) ) for i in range(7))
        print(it,'rank',rank,'localdim(incl gauge)',28-rank,'pullback-defect',f"{pull:.2e}",list(r),file=out,flush=True)
    print('nondeg',nd,file=out)
