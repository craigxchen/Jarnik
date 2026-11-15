import numpy as np, itertools, flint, sys
from scipy.optimize import linprog
P=np.array([[1,1,1,1,1,-1,-1,-1],[1,1,1,-1,-1,1,1,-1],[1,1,-1,1,-1,1,-1,1],[1,-1,1,1,-1,-1,1,1],[1,-1,-1,-1,1,1,1,1],[-1,1,1,-1,1,1,-1,1]])
D=8;L=3
d=[(P[i]-P[0])//2 for i in range(1,6)]
print("pair deltas:")
for i in range(6):
    for j in range(i+1,6):
        dl=(P[i]-P[j])//2; print(i,j,(dl>0).sum(),(dl<0).sum())
M=flint.fmpz_mat([[int(x) for x in v] for v in d]); print("rank",M.rank())
Y,nul=M.nullspace(); Y=np.array([[int(Y[r,c]) for c in range(nul)] for r in range(D)]); print("complement basis (cols)\n",Y.T)
# Gordan: exists y>0 with Y-span? need y>0 orthogonal to V, i.e. y in colspan(Y). LP: y=Y w, y>=1
res=linprog(np.zeros(nul),A_ub=-Y,b_ub=-np.ones(D),bounds=[(None,None)]*nul); print("Gordan feasible", res.status==0, res.x if res.status==0 else None)
# enumerate integer vectors in V with entries in [-4,4]
best=None; cnt=0
for a in itertools.product(range(-4,5),repeat=D):
    a=np.array(a)
    if not a.any() or (Y.T@a).any(): continue
    ap=a[a>0].sum(); am=-a[a<0].sum()
    if ap+am<=L or 2*min(ap,am)<L: print("KILL",a); break
    cnt+=1
print("in-span vectors checked",cnt)
