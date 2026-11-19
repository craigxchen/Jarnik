# Numerical probe (EVIDENCE only): Levenberg-Marquardt on the complex split system
# E(u^k - v^k)=0, k=1..5, gauge u_01=0, v_01=1; 18 complex unknowns, 20 equations.
import numpy as np, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 1)
edges=[(i,j) for i in range(5) for j in range(i+1,5)]
E=np.array([[ (1 if r in e else 0)-(1 if 4 in e else 0) for e in edges] for r in range(4)],float)
def unpack(z):
    u=np.concatenate([[0],z[:9]]); v=np.concatenate([[1],z[9:]]); return u,v
def F(z):
    u,v=unpack(z); return np.concatenate([E@(u**k-v**k) for k in range(1,6)])
def J(z):
    u,v=unpack(z); rows=[]
    for k in range(1,6):
        Ju=E*(k*u**(k-1))[None,:]; Jv=-E*(k*v**(k-1))[None,:]
        rows.append(np.hstack([Ju[:,1:],Jv[:,1:]]))
    return np.vstack(rows)
best=[]
for trial in range(int(sys.argv[2]) if len(sys.argv)>2 else 300):
    z=(rng.normal(size=18)+1j*rng.normal(size=18))
    lam=1e-2
    for it in range(400):
        f=F(z); Jz=J(z); A=Jz.conj().T@Jz; g=Jz.conj().T@f
        dz=np.linalg.lstsq(A+lam*np.eye(18),-g,rcond=None)[0]; zn=z+dz
        if not np.all(np.isfinite(zn)): break
        if np.linalg.norm(F(zn))<np.linalg.norm(f): z=zn; lam=max(lam/3,1e-12)
        else: lam*=4
        if np.linalg.norm(F(z))<1e-13: break
    if not np.all(np.isfinite(z)): continue
    u,v=unpack(z); vals=np.concatenate([u,v]); 
    d=np.abs(vals[:,None]-vals[None,:])+np.eye(20)*1e9
    scale=np.max(np.abs(vals))+1
    best.append((np.linalg.norm(F(z)), d.min()/scale))
best=np.array(best)
conv=best[best[:,0]<1e-10]
print("runs",len(best),"converged",len(conv),"min relative separation among converged:", (conv[:,1].max() if len(conv) else None))
print("residual quantiles",np.quantile(best[:,0],[0,0.1,0.5]))
