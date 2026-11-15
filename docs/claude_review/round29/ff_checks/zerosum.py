# Solve: D roots r_k in C\R such that given subsets F (list of index lists) are zero-sum: sum_{k in F} Im r_k^m = 0, m=1..L.
import numpy as np, itertools, sys
from scipy.optimize import least_squares
def resid(x, Fs, L, D):
    r=x[:D]+1j*x[D:]
    out=[]
    for F in Fs:
        rr=r[F]
        for m in range(1,L+1):
            out.append(np.sum(rr**m).imag/ (1+np.max(np.abs(rr))**m))
    # gauge: x0=0, y0=1 (real affine maps)
    out.append(10*(x[0])); out.append(10*(x[D]-1))
    return np.array(out)
def solve(Fs,L,D,tries=200,seed=0):
    rs=np.random.RandomState(seed); sols=[]
    for t in range(tries):
        x0=rs.randn(2*D)
        res=least_squares(resid,x0,args=(Fs,L,D),xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
        if res.cost<1e-24:
            r=res.x[:D]+1j*res.x[D:]
            # nondegeneracy: nonreal, distinct, no conjugate pairs
            dmin=min([abs(r[i]-r[j]) for i in range(D) for j in range(i+1,D)]+[abs(r[i]-np.conj(r[j])) for i in range(D) for j in range(D)])
            sols.append((dmin,r))
    return sols
def count_zero_sum(r,L,tol=1e-8):
    D=len(r); cnt=[]
    V=np.array([[ (rk**m).imag for rk in r] for m in range(1,L+1)])
    for mask in range(2**D):
        x=np.array([(mask>>k)&1 for k in range(D)])
        if np.max(np.abs(V@x))<tol*(1+np.abs(V).max()): cnt.append(mask)
    return cnt
if __name__=="__main__":
    D=8; L=3
    cube=list(itertools.product([0,1],repeat=3))
    Fs=[list(range(8))]+[[k for k,u in enumerate(cube) if u[i]==0] for i in range(3)]
    sols=solve(Fs,L,D,tries=300)
    print(len(sols), "solutions;", "nondegenerate:", sum(1 for d,r in sols if d>1e-3))
    sols.sort(key=lambda s:-s[0])
    for d,r in sols[:3]:
        print(d, np.round(r,6)); z=count_zero_sum(r,L); print("zero-sum subsets:",len(z))
