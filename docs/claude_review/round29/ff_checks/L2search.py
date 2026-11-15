# L=2: columns v_k=(Im rho, Im rho^2)=(y,2xy): arbitrary vectors with first coord >0; distinct roots <=> distinct columns.
# max over n distinct integer vectors of max # subsets with equal sum (sign-vector fibres = subset-sum fibres).
import itertools, numpy as np, sys
from collections import Counter
n=int(sys.argv[1]); B=int(sys.argv[2])
cand=[(a,b) for a in range(1,B+1) for b in range(-B,B+1) if np.gcd(a,abs(b))>=1]
best=0;bestv=None
masks=((np.arange(2**n)[:,None]>>np.arange(n))&1)
for comb in itertools.combinations(cand,n):
    V=np.array(comb)
    s=masks@V
    h=s[:,0]*100003+s[:,1]
    m=np.bincount(h-h.min()).max()
    if m>best: best=m;bestv=comb;print(best,comb,flush=True)
print("final",n,B,best,bestv)
