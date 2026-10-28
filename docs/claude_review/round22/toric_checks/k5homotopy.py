# Newton homotopy for the complexified k=5, d=1 system (see k5complex.py): H(z,t) = F(z) - (1-t) F(z*),
# z* random complex start. Track t: 0 -> 1 with Euler predictor + Newton corrector, adaptive step.
# Endpoints sampled; classify nondegenerate (c_j distinct nonzero, roots distinct) vs degenerate / diverging.
import numpy as np, sys
from k5complex import F, J, unpack, cvals
def track(z0, maxsteps=20000):
    b=F(z0); z=z0.copy(); t=0.0; h=0.01
    for _ in range(maxsteps):
        if t>=1: break
        h=min(h,1-t)
        # predictor: dz/dt = -J^{-1} b  (since d/dt [F(z) - (1-t)b] = J dz/dt + b = 0)
        try: dz=np.linalg.solve(J(z),-b)
        except np.linalg.LinAlgError: return z,t,'sing'
        zp=z+h*dz; tn=t+h; ok=False
        for _ in range(6):
            r=F(zp)-(1-tn)*b
            try: d=np.linalg.solve(J(zp),-r)
            except np.linalg.LinAlgError: break
            zp=zp+d
            if np.linalg.norm(d)<1e-10*(1+np.linalg.norm(zp)): ok=True; break
        if ok and np.linalg.norm(zp-z)<0.5*(1+np.linalg.norm(z)):
            z=zp; t=tn; h=min(h*1.6,0.05)
        else:
            h*=0.5
            if h<1e-12: return z,t,'stall'
        if np.linalg.norm(z)>1e6: return z,t,'diverge'
    # endgame: Newton at t=1
    for _ in range(50):
        r=F(z)
        if np.max(np.abs(r))<1e-12: break
        try: z=z+np.linalg.solve(J(z),-r)
        except np.linalg.LinAlgError: break
    return z,t,'end' if np.max(np.abs(F(z)))<1e-9 else 'endfail'
if __name__=="__main__":
    n=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); stats={}; sols=[]
    for k in range(n):
        z0=(rng.normal(size=28)+1j*rng.normal(size=28))
        z,t,st=track(z0)
        if st=='end':
            r,rp=unpack(z); c=cvals(r,rp); allr=np.concatenate([r,rp])
            dmin=min(abs(allr[a]-allr[b]) for a in range(30) for b in range(a+1,30))
            nd=min(np.min(np.abs(c)), min(abs(c[a]-c[b]) for a in range(4) for b in range(a+1,4)), dmin)
            sv=np.linalg.svd(J(z),compute_uv=False); cond=sv[-1]/sv[0]
            st='nondeg' if (nd>1e-6 and cond>1e-12) else 'deg'
            if st=='nondeg': sols.append(z); print(f"path {k}: nondeg nd={nd:.2e} cond={cond:.1e} |z|={np.linalg.norm(z):.2e}",flush=True)
        stats[st]=stats.get(st,0)+1
    print(stats)
    np.save(f"k5hom_{sys.argv[2]}.npy",np.array(sols) if sols else np.zeros((0,28)))
