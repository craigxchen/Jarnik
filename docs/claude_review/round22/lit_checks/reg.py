# reg(S_m) = max order of vanishing at e of a section of m(-K) on Lbar_k (k=5), via rank mod p.
# ord_e(sum c_u x^u) >= r  iff  c orthogonal to all polynomials of degree < r on S.
import sys, itertools, numpy as np
sys.path.insert(0,'.')
from polytope import points
def monos(n,r):
    return [e for d in range(r+1) for e in itertools.product(range(d+1),repeat=n) if sum(e)==d]
def rank_mod(A,p):
    A=A.copy()%p; nr,nc=A.shape; r=0
    for c in range(nc):
        piv=None
        nz=np.nonzero(A[r:,c])[0]
        if len(nz)==0: continue
        piv=r+nz[0]
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),p-2,p)
        A[r]=(A[r]*inv)%p
        col=A[:,c].copy(); col[r]=0
        idx=np.nonzero(col)[0]
        if len(idx): A[idx]=(A[idx]-np.outer(col[idx],A[r]))%p
        r+=1
        if r==nr: break
    return r
def reg(S,p):
    U=np.array([s[:-1] for s in S],dtype=np.int64)  # coordinates u1..u_{k-1}
    n=U.shape[1]; r=0
    while True:
        ms=monos(n,r)
        M=np.ones((len(S),len(ms)),dtype=np.int64)
        for j,e in enumerate(ms):
            col=np.ones(len(S),dtype=np.int64)
            for i,ei in enumerate(e):
                for _ in range(ei): col=(col*(U[:,i]%p))%p
            M[:,j]=col
        rk=rank_mod(M.T.copy(),p)  # rows = monomials
        if rk==len(S): return r, len(ms)
        r+=1
if __name__=='__main__':
    k=int(sys.argv[1]); p=2147483629
    for m in range(1,int(sys.argv[2])+1):
        S=points(k,m); R,nm=reg(S,p)
        print(f"k={k} m={m} |S|={len(S)} reg_p={R} ord/m={R/m:.4f} (tau bound {4-2**(3-k)})",flush=True)
