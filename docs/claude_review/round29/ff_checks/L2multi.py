# L=2 with multiplicities: multiset of n integer vectors (y>0), points = sub-multisets a in prod[0,e_k]; max fibre of sum.
import itertools, numpy as np, sys
from collections import Counter
n=int(sys.argv[1]); B=int(sys.argv[2])
cand=[(a,b) for a in range(1,B+1) for b in range(-B,B+1)]
best=0
for comb in itertools.combinations_with_replacement(range(len(cand)),n):
    cnt=Counter(comb); keys=list(cnt)
    V=np.array([cand[k] for k in keys]); E=[cnt[k] for k in keys]
    grids=np.array(list(itertools.product(*[range(e+1) for e in E])))
    s=grids@V
    h=s[:,0]*100003+s[:,1]
    m=np.bincount(h-h.min()).max()
    if m>best: best=m;print(best,[(cand[k],cnt[k]) for k in keys],flush=True)
print("final",n,B,best)
