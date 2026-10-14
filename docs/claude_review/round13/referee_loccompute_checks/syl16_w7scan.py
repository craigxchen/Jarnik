# all primitive v in {-1,0,1}^15 (first nonzero >0) with width |<v,A>| <= 7 on Syl(16,2)
import numpy as np, itertools, time
from syl16_base import COLS
A=np.load("syl16_A.npy").astype(np.float32)
supp=(A!=0).sum(1)
rng=np.random.default_rng(1)
F=np.vstack([A[supp<=4], A[rng.choice(len(A),3000,replace=False)]])
print("filter pts",len(F))
t0=time.time()
grid=np.array(list(itertools.product([-1,0,1],repeat=8)),dtype=np.float32)
grid7=np.array(list(itertools.product([-1,0,1],repeat=7)),dtype=np.float32)
surv=[]
for g in grid:
    V=np.hstack([np.repeat(g[None,:],len(grid7),0),grid7])
    nz=np.any(V!=0,axis=1)
    first=np.array([v[np.nonzero(v)[0][0]] if v.any() else 0 for v in V]) if False else None
    M=V@F.T
    ok=(M.max(1)<=3)&nz
    for v in V[ok]:
        f=v[np.nonzero(v)[0][0]]
        if f>0: surv.append(v)
print("survivors",len(surv),time.time()-t0)
S=np.array(surv,dtype=np.float32)
res=[]
for st in range(0,len(S),2000):
    M=S[st:st+2000]@A.T
    res.append(M.max(1)-M.min(1)+1)
w=np.concatenate(res)
out=S[w<=7].astype(int)
print("width<=7 dirs",len(out))
from collections import Counter
print(Counter((int(np.abs(v).sum()),int(wi)) for v,wi in zip(S.astype(int),w) if wi<=7))
np.save("syl16_w7dirs.npy",out)
# classify: columns, sigma, half-sums, quarter-sums
Cs={tuple(c) for c in COLS.astype(int)}
known=0
for v in out:
    t=tuple(int(x) for x in v)
    if t in Cs: known+=1
print("of which columns",known)
for v in out[:40]: print(v.tolist(), int(np.abs(v).sum()))
