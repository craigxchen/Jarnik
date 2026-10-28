# targeted search of the x-structured ansatz (3.4): Thm-4.1 curve in x (rho_S, S subset [3]) + conic block H = alpha + beta z
import numpy as np, itertools, sys
from scipy.optimize import least_squares
S3=[frozenset(c) for k in (1,2,3) for c in itertools.combinations((1,2,3),k)]  # 7, last = 123
MEM=[[i for i,S in enumerate(S3) if j in S] for j in (1,2,3)]
def unpack(y):
    rho=y[0:7]+1j*y[7:14]; r4=y[14]+1j*y[15]; sg=SIGN; sig=sg*np.exp(y[16]); beta=np.concatenate([[1.0],y[17:20]]); q=y[20:23]
    return rho,r4,sig,beta,q
SIGN=1.0
def F(y):
    rho,r4,sig,beta,q=unpack(y); out=[]
    for j in range(3):
        c=np.poly(rho[MEM[j]]); out+=list(c[1:4].imag)            # Thm 4.1: e1,e2,e3 real (9)
    n=np.poly(np.concatenate([rho,[r4]]))                          # monic deg 8
    out+=list(n[1:4].imag)                                          # Im n has degree <= 4 (3)
    In=n.imag; Rn=n.real
    lhs=Rn-np.polymul(In,In)[-9:]/(4*sig)+np.concatenate([np.zeros(8),[sig]])
    rhs=np.polymul(np.polymul(beta,beta),q)
    rhs=np.concatenate([np.zeros(9-len(rhs)),rhs])
    out+=list(lhs-rhs)                                              # (9)
    out+=[rho[6].real, rho[6].imag-1.0]                             # gauge: rho_123 = i
    return np.array(out)
def check(y):
    rho,r4,sig,beta,q=unpack(y)
    pts=np.concatenate([rho,[r4]]); allv=np.concatenate([pts,np.conj(pts)])
    iu=np.triu_indices(16,1); sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu])
    return sep,np.min(np.abs(pts.imag)),q[0],np.max(np.abs(pts))
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w'); nd=0; nc=0
    for it in range(N):
        SIGN=rng.choice([1.0,-1.0]); globals()['SIGN']=SIGN
        y0=np.concatenate([rng.normal(size=7)*1.5,rng.normal(size=7)*1.5,rng.normal(size=2),[rng.normal()],rng.normal(size=3),rng.normal(size=3)])
        try: sol=least_squares(F,y0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        except Exception: continue
        r=np.linalg.norm(F(sol.x))
        if not np.isfinite(r) or r>1e-10: continue
        nc+=1; sep,imin,q2,mx=check(sol.x)
        ok=sep>1e-3 and imin>1e-3 and q2>1e-6 and mx<1e4
        nd+=ok
        print(it,f"res={r:.1e} sep={sep:.1e} imin={imin:.1e} q2={q2:.2e}",'NONDEG' if ok else '',list(sol.x) if ok else '',file=out,flush=True)
    print('converged',nc,'nondeg(q2>0)',nd,file=out)
