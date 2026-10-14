# All primitive v in {-1,0,1}^15 (first nonzero > 0) with width |<v,A>| <= 7 on Syl(16,2).
# Filter on a subset of A (a necessary condition), then exact check on all of A (float32 exact for these ints).
import numpy as np, itertools, time
from syl16_base import COLS
A=np.load("syl16_A.npy").astype(np.float32)
supp=(A!=0).sum(1)
rng=np.random.default_rng(1)
F=np.vstack([A[supp<=3], A[rng.choice(len(A),1500,replace=False)]])
print("filter pts",len(F),flush=True)
t0=time.time()
g10=np.array(list(itertools.product([-1,0,1],repeat=10)),dtype=np.float32)
g5=np.array(list(itertools.product([-1,0,1],repeat=5)),dtype=np.float32)
FT1=F[:,:5].T; FT2=F[:,5:].T
P2=g10@FT2
surv=[]
for g in g5:
    p1=g@FT1
    M=P2+p1[None,:]
    ok=np.where(M.max(1)<=3)[0]
    for i in ok:
        v=np.concatenate([g,g10[i]])
        if v.any() and v[np.nonzero(v)[0][0]]>0: surv.append(v)
print("survivors",len(surv),time.time()-t0,flush=True)
S=np.array(surv,dtype=np.float32)
ws=[]
for st in range(0,len(S),500):
    M=S[st:st+500]@A.T
    ws.append(M.max(1)-M.min(1)+1)
w=np.concatenate(ws)
out=S[w<=7].astype(int); wo=w[w<=7]
print("width<=7 dirs",len(out),time.time()-t0)
from collections import Counter
print(Counter((int(np.abs(v).sum()),int(x)) for v,x in zip(out,wo)))
np.save("syl16_w7dirs.npy",out); np.save("syl16_w7w.npy",wo)
