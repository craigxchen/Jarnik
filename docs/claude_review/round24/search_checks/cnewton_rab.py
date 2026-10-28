# complexified system + Rabinowitsch w_j*kappa_j = 1 ; Newton from random starts
import numpy as np, sys
from cnewton import *
def FJ(Z,L,b):
    z=Z[:30]; w=Z[30:]
    F0,J0=F_J(z,L,b)   # 30 x 30 (28 eqs + 2 slices)
    F=np.zeros(34,dtype=complex); J=np.zeros((34,34),dtype=complex)
    F[:30]=F0; J[:30,:30]=J0
    r=z[:15]; rp=z[15:]
    for j in range(4):
        c,Jc=coef_and_jac(r[MEM[j]]); cp,Jcp=coef_and_jac(rp[MEM[j]])
        kj=c[4]-cp[4]
        F[30+j]=w[j]*kj-1
        J[30+j,MEM[j]]=w[j]*Jc[4]; J[30+j,[15+m for m in MEM[j]]]=-w[j]*Jcp[4]; J[30+j,30+j]=kj
    return F,J
if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); out=open(sys.argv[3],'w'); rng=np.random.default_rng(seed)
    nc=nd=0
    for it in range(N):
        L=rng.normal(size=(2,30))+1j*rng.normal(size=(2,30))
        z=(rng.normal(size=30)+1j*rng.normal(size=30))*rng.choice([0.5,1,2])
        Z=np.concatenate([z,1/kap(z)])
        b=L@z
        conv=False
        for step in range(400):
            F,J=FJ(Z,L,b)
            try: dZ=np.linalg.solve(J,-F)
            except np.linalg.LinAlgError: break
            lam=1.0; n0=np.linalg.norm(F)
            while lam>1e-5:
                Zn=Z+lam*dZ
                if np.linalg.norm(FJ(Zn,L,b)[0])<n0*(1-0.3*lam): break
                lam/=2
            Z=Zn
            if np.linalg.norm(FJ(Z,L,b)[0])<1e-11: conv=True; break
            if np.max(np.abs(Z))>1e8: break
        if not conv: continue
        nc+=1; z=Z[:30]
        iu=np.triu_indices(30,1); sep=np.min(np.abs(z[:,None]-z[None,:])[iu])
        k=np.min(np.abs(kap(z))); mn=np.min(np.abs(z)); mx=np.max(np.abs(z))
        ok= sep>1e-4 and k>1e-6 and mn>1e-4 and mx<1e4
        nd+=ok
        print(it,f"sep={sep:.2e} kap={k:.2e} min={mn:.2e} max={mx:.2e}",'NONDEG' if ok else '', list(z) if ok else '',file=out,flush=True)
    print('converged',nc,'nondeg',nd,file=out)
