# count zero-sum sets at a refined solution, 60 digits; prints count and the smallest nonzero residual
import sys, itertools, numpy as np, mpmath as mp
sys.path.insert(0,'.'); from certify import certify
mp.mp.dps=60
def count(z,D,L):
    q=[mp.mpc(z[i],z[D+i]) for i in range(D)]
    P=[[mp.im(qk**m) for m in range(1,L+1)] for qk in q]
    zs=[];mn=mp.mpf(10)
    for F in itertools.product([0,1],repeat=D):
        v=[sum(P[k][m] for k in range(D) if F[k]) for m in range(L)]
        nv=mp.sqrt(sum(x*x for x in v))
        if nv<mp.mpf(10)**-40: zs.append(F)
        elif nv<mn: mn=nv
    return zs,mn
if __name__=="__main__":
    import itertools
    D=12;L=5;U=list(itertools.product([0,1],repeat=3)); cols=[]
    for u in U: cols+=[u]*((12+4*(-1)**sum(u))//8)
    cols=np.array(cols); Fs=[[1]*D]+[list(1-cols[:,i]) for i in range(3)]
    qs=np.load("hc12_lam4.npy")
    ok,res,qq,fix,z=certify(qs[0],Fs,L)
    zs,mn=count(z,D,L); print("certified",ok,"zero-sum sets",len(zs),"rank",np.linalg.matrix_rank(np.array(zs)),"min nonzero residual",mp.nstr(mn,5))
    print("q =",[ (mp.nstr(z[i],20),mp.nstr(z[D+i],20)) for i in range(D)])
    print("fixed",[(j,str(v)) for j,v in fix])
