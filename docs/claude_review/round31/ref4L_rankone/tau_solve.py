# tau_e = omega * z^{B_e} (all vertex products 1); equations sum_{e at i}(tau_e+beta/tau_e)=0
import numpy as np, sys, mpmath as mp
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
B=np.array([[0,1,-1,0,-1,1,0,0,0,0],[0,1,0,-1,-1,0,1,0,0,0],[1,0,-1,0,-1,0,0,1,0,0],[1,0,0,-1,-1,0,0,0,1,0],[1,1,-1,-1,-1,0,0,0,0,1]]).T  # 10x5
rng=np.random.default_rng(int(sys.argv[1]))
def tau(lz,om): return om*np.exp(B@lz)
def F(lz,beta,om):
    t=tau(lz,om); return np.array([np.sum(t[inc[i]]+beta/t[inc[i]]) for i in range(5)])
def Jac(lz,beta,om):
    t=tau(lz,om); D=(t-beta/t)[:,None]*B  # d tau/d lz = tau*B
    return np.array([D[inc[i]].sum(0) for i in range(5)])
res={}
for beta in [0.7+0.3j, -1.3+0.2j, 2.1-1.1j]:
  for om in [1.0, 1j]:
    sols=[]
    for t in range(4000):
        lz=rng.normal(size=5)+1j*rng.normal(size=5)*np.pi
        for it in range(80):
            f=F(lz,beta,om)
            if np.max(np.abs(f))<1e-13: break
            try: lz=lz-np.linalg.solve(Jac(lz,beta,om),f)
            except np.linalg.LinAlgError: break
            if np.max(np.abs(lz.real))>30: break
        if np.max(np.abs(F(lz,beta,om)))<1e-11:
            tt=tau(lz,om)
            key=tuple(np.round(np.sort_complex(np.round(tt,8)),6))
            if not any(np.allclose(np.sort_complex(tt),np.sort_complex(np.array(s)),atol=1e-6) for s in sols):
                sols.append(tt)
    nd=0
    for tt in sols:
        rho=tt-beta/tt
        g=min(abs(rho[a]-rho[b]) for a in range(10) for b in range(a))
        if g>1e-6: nd+=1
    print('beta',beta,'omega',om,': distinct tau-solutions',len(sols),' nondegenerate (all 10 rho distinct):',nd,flush=True)
