import numpy as np, itertools, sys
sys.path.insert(0,'.'); from basesolve import *; from certify import *
D=12; L=5
U=list(itertools.product([0,1],repeat=3)); lam=4
cols=[]
for u in U: cols+=[u]*((12+lam*(-1)**sum(u))//8)
cols=np.array(cols)
Fs=[[1]*D]+[list(1-cols[:,i]) for i in range(3)]
qs=np.load("hc12_lam4.npy")
for q in qs[:2]:
    ok,res,qq,fix,z=certify(q,Fs,L)
    qq=np.array(qq); print("Krawczyk ok:",ok,"residual %.1e"%float(res),"fixed",[(j,str(v)) for j,v in fix])
    print("min |Im q| %.4f  min root separation %.4f"%(np.abs(qq.imag).min(),degen(qq)))
    S=zero_sum_sets(qq,L,1e-20); print("zero-sum sets at 60 digits (tol 1e-20, double eval):",len(S))
