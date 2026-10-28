# k=5, block degree d=1 balanced extremal curves  <=>  for j=1..4 the 8-multisets R_j={r_T: T∋j} and conj(R_j)
# have equal power sums p_1..p_7 (ideal complex PTE), i.e. Im sum_{T∋j} r_T^n = 0, n=1..7 (28 eqs).
# Nondegeneracy (Rabinowitsch): y_T*Im r_T = 1 (15), w_j*c_j = 1 (4), c_j = Im prod_{T∋j}(-r_T).
# Gauge (real affine s -> lam s + mu): sum Re r_T = 0, Im r_top = 1.  Square: 47 unknowns, 47 eqs. Analytic Jacobian.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
Ts=[T for r in range(1,5) for T in itertools.combinations(range(4),r)]
iTop=Ts.index((0,1,2,3)); idx=[i for i in range(15) if i!=iTop]
Mem=np.array([[1.0 if j in T else 0.0 for T in Ts] for j in range(4)])   # 4 x 15
def unpack(z):
    a=np.empty(15); b=np.empty(15)
    a[:14]=z[:14]; a[14]=-z[:14].sum(); b[idx]=z[14:28]; b[iTop]=1.0
    return a+1j*b, z[28:43], z[43:47]
# d(a)/dz, d(b)/dz
Da=np.zeros((15,47)); Da[:14,:14]=np.eye(14); Da[14,:14]=-1
Db=np.zeros((15,47)); 
for k,i in enumerate(idx): Db[i,14+k]=1
def F(z):
    r,y,w=unpack(z)
    pw=np.array([r**n for n in range(1,8)])          # 7 x 15
    E1=(Mem@pw.imag.T).T.ravel()                      # 28
    E2=y*r.imag-1
    P0=np.array([np.prod(-r[Mem[j]>0]) for j in range(4)]); c=P0.imag
    E3=w*c-1
    return np.concatenate([E1,E2,E3])
def J(z):
    r,y,w=unpack(z); out=np.zeros((47,47))
    dr=np.array([n*r**(n-1) for n in range(1,8)])    # d r^n / d r (holomorphic), 7x15
    # Im f: d/da = Im f', d/db = Re f'
    dIa=dr.imag; dIb=dr.real
    row=0
    for n in range(7):
        for j in range(4):
            out[row]=(Mem[j]*dIa[n])@Da+(Mem[j]*dIb[n])@Db; row+=1
    for t in range(15):
        out[row]=y[t]*Db[t]; out[row,28+t]+=r[t].imag; row+=1
    for j in range(4):
        sel=Mem[j]>0; P0=np.prod(-r[sel]); c=P0.imag
        g=np.where(sel,P0/r,0)        # d P0/d r_T = P0/r_T for T∋j
        out[row]=w[j]*(g.imag@Da+g.real@Db); out[row,43+j]=c; row+=1
    return out
def run(ntry,seed,verbose=True):
    rng=np.random.default_rng(seed); sols=[]
    for it in range(ntry):
        sc=rng.choice([0.3,1,3])
        a=rng.normal(size=14)*sc; b=rng.normal(size=14)*sc
        r,_,_=unpack(np.concatenate([a,b,np.ones(19)]))
        y=1/r.imag; c=np.array([np.prod(-r[Mem[j]>0]).imag for j in range(4)])
        z0=np.concatenate([a,b,y,1/c])
        sol=least_squares(F,z0,jac=J,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        m=np.max(np.abs(sol.fun))
        if m<1e-10:
            r,y,w=unpack(sol.x)
            sv=np.linalg.svd(J(sol.x),compute_uv=False)
            sols.append(sol.x)
            if verbose: print(f"try {it}: res {m:.1e} min|Im r|={np.min(abs(r.imag)):.3e} c={np.round(1/w,5)} sv_min/max={sv[-1]/sv[0]:.2e}",flush=True)
    return sols
if __name__=="__main__":
    n,seed=int(sys.argv[1]),int(sys.argv[2])
    s=run(n,seed); print("found",len(s),"of",n)
    np.save(f"k5pte_{seed}.npy",np.array(s) if s else np.zeros((0,47)))
