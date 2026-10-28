# complexified system: rho_T, rho'_T independent; e_k(R_j) = e_k(R'_j), k != 4; random complex slices; Newton
import itertools, numpy as np, sys
SUBS=[frozenset(c) for k in range(1,5) for c in itertools.combinations(range(1,5),k)]
IDX={T:i for i,T in enumerate(SUBS)}; MEM=[[IDX[T] for T in SUBS if j in T] for j in range(1,5)]
KS=[1,2,3,5,6,7,8]
def coef_and_jac(roots):
    c=np.poly(roots)  # length 9
    J=np.zeros((9,len(roots)),dtype=complex)
    for i in range(len(roots)):
        other=np.poly(np.delete(roots,i))  # deg 7
        J[1:,i]=-other
    return c,J
def F_J(z,L,b):
    r=z[:15]; rp=z[15:]
    F=[]; J=np.zeros((30,30),dtype=complex); row=0
    for j in range(4):
        c,Jc=coef_and_jac(r[MEM[j]]); cp,Jcp=coef_and_jac(rp[MEM[j]])
        for k in KS:
            F.append(c[k]-cp[k])
            J[row,[m for m in MEM[j]]]=Jc[k]; J[row,[15+m for m in MEM[j]]]=-Jcp[k]; row+=1
    F=np.array(F+list(L@z-b)); J[28:]=L
    return F,J
def kap(z):
    r=z[:15]; rp=z[15:]
    return np.array([np.poly(r[MEM[j]])[4]-np.poly(rp[MEM[j]])[4] for j in range(4)])
if __name__=='__main__':
    seed=int(sys.argv[1]); N=int(sys.argv[2]); out=open(sys.argv[3],'w'); rng=np.random.default_rng(seed)
    nc=nd=0
    for it in range(N):
        L=rng.normal(size=(2,30))+1j*rng.normal(size=(2,30))
        z=(rng.normal(size=30)+1j*rng.normal(size=30))
        b=L@z + 0*(rng.normal(size=2))
        conv=False
        for step in range(300):
            F,J=F_J(z,L,b)
            try: dz=np.linalg.solve(J,-F)
            except np.linalg.LinAlgError: break
            lam=1.0
            n0=np.linalg.norm(F)
            while lam>1e-4:
                zn=z+lam*dz
                if np.linalg.norm(F_J(zn,L,b)[0])<n0*(1-0.3*lam): break
                lam/=2
            z=zn
            if np.linalg.norm(F_J(z,L,b)[0])<1e-11: conv=True; break
        if not conv: continue
        nc+=1
        allv=z; iu=np.triu_indices(30,1); sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu])
        k=np.min(np.abs(kap(z))); mn=np.min(np.abs(z)); mx=np.max(np.abs(z))
        sv=np.linalg.svd(F_J(z,L,b)[1],compute_uv=False); cond=sv[-1]/sv[0]
        ok= sep>1e-4 and k>1e-6 and mn>1e-4 and mx<1e4
        nd+=ok
        print(it,f"sep={sep:.2e} kap={k:.2e} min={mn:.2e} max={mx:.2e} smin/smax={cond:.2e}",'NONDEG' if ok else '', list(z) if ok else '',file=out,flush=True)
    print('converged',nc,'nondeg',nd,file=out)
