# Search for balanced Liouville-extremal rational curves at k=5 of minimal degree (Gaussian linear blocks, d=1):
# roots r_T (T nonempty subset of [4]), P_j(s) = prod_{T ∋ j} (s - r_T), j=1..4 (deg 8);
# need Im(coeff of s^i) = 0 for i=1..7 (28 real eqs). Gauge: sum Re r_T = 0, Im r_{1234} = 1.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
Ts=[T for r in range(1,5) for T in itertools.combinations(range(4),r)]
iTop=Ts.index((0,1,2,3))
def unpack(z):
    a=np.zeros(15); b=np.zeros(15)
    a[:14]=z[:14]; a[14]=-z[:14].sum()
    idx=[i for i in range(15) if i!=iTop]
    b[idx]=z[14:]; b[iTop]=1.0
    return a+1j*b
def polys(r):
    out=[]
    for j in range(4):
        P=np.array([1.0+0j])
        for t,T in enumerate(Ts):
            if j in T: P=np.convolve(P,np.array([1.0,-r[t]]))
        out.append(P)   # highest degree first, length 9
    return out
def F(z):
    r=unpack(z); res=[]
    for P in polys(r): res.extend(P[1:8].imag)   # coefficients of s^7..s^1
    return np.array(res)
def run(ntry,seed):
    rng=np.random.default_rng(seed); sols=[]
    for it in range(ntry):
        z0=rng.normal(size=28)*rng.choice([0.5,1,2,4])
        sol=least_squares(F,z0,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=20000)
        if np.max(np.abs(sol.fun))<1e-10:
            r=unpack(sol.x); Ps=polys(r)
            c=[P[8].imag for P in Ps]
            J=sol.jac; sv=np.linalg.svd(J,compute_uv=False)
            sols.append((min(abs(r.imag)),c,sv[-1]/sv[0],r))
            print(f"try {it}: res {np.max(np.abs(sol.fun)):.1e} min|Im r|={min(abs(r.imag)):.3e} c={np.round(c,6)} cond^-1={sv[-1]/sv[0]:.2e}",flush=True)
    return sols
if __name__=="__main__":
    sols=run(int(sys.argv[1]),int(sys.argv[2]))
    np.save(f"k5sols_{sys.argv[2]}.npy",np.array([s[3] for s in sols]))
