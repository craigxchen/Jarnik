import numpy as np, itertools, time, sys
from flint import nmod_mat
from syl16_base import COLS, pts, build
from reng import REng, PR
from rlib import J_exact, monos
A=build(1).astype(np.int64); print("|A|",len(A))
d=15; DT=4
JA=J_exact([tuple(int(x) for x in p) for p in A],DT,cap=2)
print("dim J(A)_D exact",[len(JA[D]) for D in range(DT+1)],flush=True)
g10=np.array(list(itertools.product([-1,0,1],repeat=10)),dtype=np.float32)
g5=np.array(list(itertools.product([-1,0,1],repeat=5)),dtype=np.float32)
Af=A.astype(np.float32)
P2=g10@Af[:,5:].T
dirs=[]
for g in g5:
    M=P2+(g@Af[:,:5].T)[None,:]
    w=M.max(1)-M.min(1)+1
    for i in np.where(w<=5)[0]:
        v=np.concatenate([g,g10[i]]).astype(int)
        if v.any() and v[np.nonzero(v)[0][0]]>0: dirs.append(tuple(v))
print("dirs width<=5:",len(dirs),flush=True)
from collections import Counter
print(Counter((sum(abs(x) for x in v),len(np.unique(A@np.array(v)))) for v in dirs))
E=REng(d,DT)
for DTT in [3,4]:
    E=REng(d,DTT)
    rows=[]
    for v in dirs:
        s=len(np.unique(A@np.array(v)))
        if s<=DTT:
            g=E.pullback([1],[(s,)],[list(v)])
            if g.any(): rows+=E.ideal_rows(g,s,DTT)
    M=np.array(rows); r=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR).rank()
    print("D",DTT,"J_w(width<=5 dirs)",r,"J(A)",len(JA[DTT]),flush=True)
    t0=time.time(); npair=0
    for U in itertools.combinations(dirs,2):
        Ua=np.array(U); img=np.unique(A@Ua.T,axis=0)
        if len(np.unique(img[:,0]))*len(np.unique(img[:,1]))==len(img): continue
        npair+=1
        J=J_exact([tuple(int(x) for x in p) for p in img],DTT)
        for D in range(DTT+1):
            for vec in J[D]:
                g=E.pullback(vec,monos(2,D),[list(U[0]),list(U[1])])
                if g.any(): rows+=E.ideal_rows(g,D,DTT)
    M=np.array(rows); r=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR).rank()
    print("D",DTT,"+ all non-grid pairs of these dirs (%d)"%npair,r,"J(A)",len(JA[DTT]),time.time()-t0,flush=True)
