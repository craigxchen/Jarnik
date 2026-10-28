# calibration on the SPLIT analogue (split Thm 2 class): 15 real places rho_T (T subset [4]) with orientation o_T;
# u_j = N_j/D_j, N_j = prod_{T∋j, o_T=+}(s-rho_T), D_j = prod_{T∋j, o_T=-}(s-rho_T) (4 and 4); contact (2,2): N_j - D_j = kappa s^2
import numpy as np, itertools, sys
from scipy.optimize import least_squares
exec(open('signpat.py').read().split('print(len(pats))')[0])
SUBS_=SUBS; IDX={T:i for i,T in enumerate(SUBS_)}
MEM=[[IDX[T] for T in SUBS_ if j in T] for j in range(1,5)]
def F(x,o,L,c0):
    out=[]
    for j in range(4):
        num=[x[m] for m in MEM[j] if o[m]>0]; den=[x[m] for m in MEM[j] if o[m]<0]
        d=np.poly(num)-np.poly(den)   # coefficients s^4..s^0, s^4 cancels
        out+= [d[1],d[3],d[4]]
    out.append(np.log(abs(x[14])+1e-300))
    return np.concatenate([out,0.1*(L@x-c0)])
def nondeg(x):
    iu=np.triu_indices(15,1); return np.min(np.abs(x[:,None]-x[None,:])[iu]), np.min(np.abs(x)), np.max(np.abs(x))
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); nc=nd=0
    for it in range(N):
        o=np.array(pats[rng.integers(len(pats))]); L=rng.normal(size=(2,15)); c0=rng.normal(size=2)
        x0=rng.normal(size=15)*rng.choice([0.5,1,2])
        try: sol=least_squares(F,x0,args=(o,L,c0),method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        except Exception: continue
        r=np.linalg.norm(F(sol.x,o,L,c0)[:12])
        if r>1e-10: continue
        nc+=1; sep,mn,mx=nondeg(sol.x)
        # also kappa != 0
        ok=sep>1e-3 and mn>1e-4 and mx<1e4
        nd+=ok
    print('split (2,2): starts',N,'converged',nc,'nondegenerate',nd)
