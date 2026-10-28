# test whether Poly_{<=r} spans R^{S_m} (k=5) mod p; if yes, every section of m(-K) has ord_e <= r.
import sys, numpy as np, time
sys.path.insert(0,'.')
from polytope import points
from reg import monos, rank_mod
k,m,r=map(int,sys.argv[1:4]); p=2147483629
t=time.time(); S=points(k,m); print('|S|',len(S),time.time()-t,flush=True)
U=np.array([s[:-1] for s in S],dtype=np.int64)%p
ms=monos(U.shape[1],r)
M=np.ones((len(S),len(ms)),dtype=np.int64)
for j,e in enumerate(ms):
    col=np.ones(len(S),dtype=np.int64)
    for i,ei in enumerate(e):
        for _ in range(ei): col=(col*U[:,i])%p
    M[:,j]=col
rk=rank_mod(M,p)
print(f"k={k} m={m} r={r} |S|={len(S)} #monos={len(ms)} rank_p={rk} spans={rk==len(S)} => ord/m <= {r/m if rk==len(S) else 'n/a'}",flush=True)
