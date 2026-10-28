# fast LM with analytic Jacobian; unknowns: rho_T = a_T + i*eps_T*exp(b_T) (sign pattern) ; equations Im e_k(R_j), k in K
import itertools, numpy as np, sys
from scipy.optimize import least_squares
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
IDX={T:i for i,T in enumerate(SUBS)}; MEM=[np.array([IDX[T] for T in SUBS if j in T]) for j in range(1,5)]
def polyjac(roots):
    c=np.poly(roots); n=len(roots); J=np.zeros((n+1,n),dtype=complex)
    for i in range(n):
        J[1:,i]=-np.poly(np.delete(roots,i))
    return c,J
class Prob:
    def __init__(s,K,eps,nslice,rng):
        s.K=np.array(K); s.eps=np.array(eps,dtype=float); s.L=rng.normal(size=(nslice,30)); s.c0=rng.normal(size=nslice)
    def rho(s,x): return x[:15]+1j*s.eps*np.exp(x[15:])
    def F(s,x):
        r=s.rho(x); out=[]
        for j in range(4):
            c=np.poly(r[MEM[j]]); out.append(c[s.K].imag)
        g=np.log(np.abs(r[14]))
        return np.concatenate(out+[[g],0.1*(s.L@x-s.c0)])
    def J(s,x):
        r=s.rho(x); rows=[]
        dr_da=np.ones(15,dtype=complex); dr_db=1j*s.eps*np.exp(x[15:])
        for j in range(4):
            c,Jc=polyjac(r[MEM[j]])
            blk=np.zeros((len(s.K),30))
            sub=Jc[s.K]  # d c_k / d rho_m (holomorphic)
            blk[:,MEM[j]]=(sub*dr_da[MEM[j]]).imag
            blk[:,15+MEM[j]]=(sub*dr_db[MEM[j]]).imag
            rows.append(blk)
        g=np.zeros((1,30)); g[0,14]=(np.conj(r[14])*dr_da[14]).real/abs(r[14])**2; g[0,29]=(np.conj(r[14])*dr_db[14]).real/abs(r[14])**2
        return np.vstack(rows+[g,0.1*s.L])
    def eqnorm(s,x):
        return np.linalg.norm(s.F(x)[:4*len(s.K)])
def diag(r):
    allv=np.concatenate([r,np.conj(r)]); iu=np.triu_indices(30,1)
    D=np.abs(allv[:,None]-allv[None,:])[iu]
    return dict(imin=np.min(np.abs(r.imag)),sep=np.min(D),rmax=np.max(np.abs(r)),rmin=np.min(np.abs(r)))
def isgood(r,res):
    d=diag(r); return res<1e-10 and d['imin']>1e-3 and d['sep']>1e-3 and d['rmax']<1e4 and d['rmin']>1e-4
def Kof(mi,m0): return list(range(1,mi))+list(range(9-m0,9))
def run_direct(K,eps,rng,gauge,x0=None,nfev=3000):
    ns=max(1,30-1-4*len(K)) if gauge==2 else max(1,30-1-4*len(K))
    P=Prob(K,eps,ns,rng)
    if x0 is None: x0=np.concatenate([rng.normal(size=15)*rng.choice([0.5,1,2]), rng.normal(size=15)*0.7])
    sol=least_squares(P.F,x0,jac=P.J,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=nfev)
    return sol.x,P.eqnorm(sol.x),P
