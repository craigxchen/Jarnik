# complex damped Newton search for nondegenerate solutions of square type systems (root formulation r_k,s_k)
import numpy as np, pickle, sys, itertools
from cuttypes import CUTS
from gb import cuts_to_columns
def make(cols,h,npts=5,slices=0,seed=0):
    n=len(cols)
    C=np.array([[a[i]-a[0] for (a,e) in cols] for i in range(1,npts)],dtype=float)  # (npts-1) x n
    E=np.array([e for a,e in cols],dtype=float)
    rng=np.random.default_rng(seed+999)
    S=rng.normal(size=(slices,2*n))+1j*rng.normal(size=(slices,2*n)); s0=rng.normal(size=slices)+1j*rng.normal(size=slices)
    def F(z):
        r,s=z[:n],z[n:]
        out=[]
        for m in range(1,h):
            out.append(C@(r**m-s**m))
        out.append([np.sum(E*(r+s))]); out.append([r[0]-s[0]-1])
        if slices: out.append(S@z+s0)
        return np.concatenate(out)
    def J(z):
        r,s=z[:n],z[n:]
        rows=[]
        for m in range(1,h):
            rows.append(np.hstack([C*(m*r**(m-1)),-C*(m*s**(m-1))]))
        g=np.concatenate([E,E])[None,:]; rows.append(g)
        g2=np.zeros((1,2*n),complex); g2[0,0]=1; g2[0,n]=-1; rows.append(g2)
        if slices: rows.append(S)
        return np.vstack(rows)
    return F,J,n
def search(cols,h,starts=200,npts=5,slices=0,seed=0,verbose=False):
    F,J,n=make(cols,h,npts,slices,seed)
    rng=np.random.default_rng(seed)
    found=[]; degenerate=0; fail=0
    for st in range(starts):
        z=(rng.normal(size=2*n)+1j*rng.normal(size=2*n))*rng.uniform(0.3,3)
        ok=False
        for it in range(200):
            f=F(z)
            nf=np.linalg.norm(f)
            if nf<1e-13: ok=True;break
            try:
                dz=np.linalg.lstsq(J(z),-f,rcond=None)[0]
            except Exception: break
            lam=1.0
            while lam>1e-6:
                zn=z+lam*dz
                if np.linalg.norm(F(zn))<nf: break
                lam/=2
            z=zn
            if not np.all(np.isfinite(z)) or np.max(np.abs(z))>1e8: break
        if not ok: fail+=1; continue
        roots=z; sc=max(1,np.max(np.abs(roots)))
        dmin=min(abs(a-b) for a,b in itertools.combinations(roots,2))/sc
        if dmin<1e-5: degenerate+=1
        else: found.append((dmin,z))
    return found,degenerate,fail
if __name__=="__main__":
    D,h=int(sys.argv[1]),int(sys.argv[2]); starts=int(sys.argv[3]); sl=int(sys.argv[4]) if len(sys.argv)>4 else 0
    js=sys.argv[5] if len(sys.argv)>5 else "all"
    T=pickle.load(open(f'types_{D}_{h}.pkl','rb'))
    rng_=range(len(T)) if js=="all" else [int(x) for x in js.split(",")]
    for j in rng_:
        cols=cuts_to_columns(T[j],CUTS)
        found,deg,fail=search(cols,h,starts,slices=sl,seed=j)
        print(D,h,j,"nondeg",len(found),"degenerate",deg,"nonconv",fail, "best sep",max([f[0] for f in found],default=None),flush=True)
