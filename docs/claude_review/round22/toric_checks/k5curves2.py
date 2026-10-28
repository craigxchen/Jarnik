# As k5curves.py, with Rabinowitsch variables excluding degenerate components:
# y_T * Im r_T = 1 (15 eqs), w_j * c_j = 1 (4 eqs), c_j = Im P_j(0). Square system 47 x 47.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
Ts=[T for r in range(1,5) for T in itertools.combinations(range(4),r)]
iTop=Ts.index((0,1,2,3)); idx=[i for i in range(15) if i!=iTop]
def unpack(z):
    a=np.zeros(15); b=np.zeros(15)
    a[:14]=z[:14]; a[14]=-z[:14].sum()
    b[idx]=z[14:28]; b[iTop]=1.0
    return a+1j*b, z[28:43], z[43:47]
def polys(r):
    out=[]
    for j in range(4):
        P=np.array([1.0+0j])
        for t,T in enumerate(Ts):
            if j in T: P=np.convolve(P,np.array([1.0,-r[t]]))
        out.append(P)
    return out
def F(z):
    r,y,w=unpack(z); res=[]; Ps=polys(r)
    for P in Ps: res.extend(P[1:8].imag)
    res.extend(y*r.imag-1)
    res.extend(w*np.array([P[8].imag for P in Ps])-1)
    return np.array(res)
def run(ntry,seed):
    rng=np.random.default_rng(seed); found=[]
    best=1e9
    for it in range(ntry):
        a=rng.normal(size=14)*rng.choice([0.5,1,2,4]); b=rng.normal(size=14)*rng.choice([0.5,1,2])
        r,_,_=unpack(np.concatenate([a,b,np.ones(19)]))
        y=1/r.imag; c=np.array([P[8].imag for P in polys(r)]); w=1/c
        z0=np.concatenate([a,b,y,w])
        sol=least_squares(F,z0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=40000)
        m=np.max(np.abs(sol.fun)); best=min(best,m)
        if m<1e-9:
            r,y,w=unpack(sol.x); c=1/w
            print(f"try {it}: res {m:.1e} min|Im r|={min(abs(r.imag)):.3e} c={np.round(c,6)}",flush=True)
            found.append(sol.x)
    print("best residual",best)
    return found
if __name__=="__main__":
    f=run(int(sys.argv[1]),int(sys.argv[2]))
    np.save(f"k5sols2_{sys.argv[2]}.npy",np.array(f))
