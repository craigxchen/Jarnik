import numpy as np, itertools, sys
sys.path.insert(0,'.'); from basesolve import *
D=12; L=5
U=list(itertools.product([0,1],repeat=3))
for lam in (4,-4):
    cols=[]
    for u in U:
        m=(12+lam*(-1)**sum(u))//8; cols+= [u]*m
    cols=np.array(cols)  # D x 3
    Fs=[np.ones(D,int)]+[1-cols[:,i] for i in range(3)]+[cols[:,i] for i in range(3)]
    Fs=np.array(Fs)
    out=solve_base(Fs,L,starts=300,seed=1)
    good=[o for o in out if o[0]<1e-10 and o[1]>1e-3]
    print("lambda",lam,"converged nondegenerate:",len(good),"best cost %.2e degen %.3e"%(out[0][0],out[0][1]))
    if good:
        q=good[0][2]; S=zero_sum_sets(q,L); 
        print("zero-sum sets:",len(S),"rank",np.linalg.matrix_rank(S)); print(np.round(q,6))
        np.save(f"hc12_lam{lam}.npy",np.array([g[2] for g in good]))
