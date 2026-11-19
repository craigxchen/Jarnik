# real Newton search: x_k,Y_k real, rho_k=x_k+iY_k, equations sum_k c_k Im(rho_k^m)=0, gauge sum x=0, Y_0=1
import numpy as np, pickle, sys, itertools
from cuttypes import CUTS
from gb import cuts_to_columns
def search(cols,h,starts,seed=0,npts=5):
    n=len(cols)
    C=np.array([[a[i]-a[0] for (a,e) in cols] for i in range(1,npts)],dtype=float)
    E=np.array([e for a,e in cols],float)
    def F(v):
        x,Y=v[:n],v[n:]; rho=x+1j*Y
        out=[C@np.imag(rho**m) for m in range(1,h)]
        out.append([E@x]); out.append([Y[0]-1])
        return np.concatenate(out)
    def J(v):
        x,Y=v[:n],v[n:]; rho=x+1j*Y
        rows=[]
        for m in range(1,h):
            d=m*rho**(m-1)  # d rho^m/dx = m rho^{m-1}; d/dY = i m rho^{m-1}
            rows.append(np.hstack([C*np.imag(d),C*np.imag(1j*d)]))
        g=np.zeros((1,2*n)); g[0,:n]=E; rows.append(g)
        g=np.zeros((1,2*n)); g[0,n]=1; rows.append(g)
        return np.vstack(rows)
    rng=np.random.default_rng(seed); found=[]; deg=0; fail=0
    for st in range(starts):
        v=rng.normal(size=2*n)*rng.uniform(0.3,3)
        ok=False
        for it in range(300):
            f=F(v); nf=np.linalg.norm(f)
            if nf<1e-13: ok=True;break
            dz=np.linalg.lstsq(J(v),-f,rcond=None)[0]
            lam=1.0
            while lam>1e-7:
                vn=v+lam*dz
                if np.linalg.norm(F(vn))<nf: break
                lam/=2
            v=vn
            if not np.all(np.isfinite(v)) or np.max(np.abs(v))>1e8: break
        if not ok: fail+=1; continue
        x,Y=v[:n],v[n:]; rho=x+1j*Y; roots=np.concatenate([rho,np.conj(rho)])
        sc=max(1,np.max(np.abs(roots)))
        sep=min(abs(a-b) for a,b in itertools.combinations(roots,2))/sc
        if sep<1e-5: deg+=1
        else: found.append((sep,v))
    return found,deg,fail
if __name__=="__main__":
    D,h,starts=map(int,sys.argv[1:4]); js=sys.argv[4] if len(sys.argv)>4 else "all"
    T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
    rng_=range(len(T)) if js=="all" else [int(x) for x in js.split(",")]
    for j in rng_:
        cols=cuts_to_columns(T[j],CUTS)
        found,deg,fail=search(cols,h,starts,seed=j)
        print(D,h,j,"REAL nondeg",len(found),"degenerate",deg,"nonconv",fail,"best sep",max([f[0] for f in found],default=None),flush=True)
        if found: pickle.dump((cols,found),open(f'real_{D}_{h}_{j}.pkl','wb'))
