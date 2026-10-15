# saturated rank of each feasible 7-multiset in r_rest output; pools given in order
import numpy as np, pickle, sys
from flint import nmod_mat
from collections import Counter
d=pickle.load(open('r_box.pkl','rb')); D=d['D'].astype(np.int64); tau=d['tau']
PR=2147483647
E=[(a,b) for a in range(12) for b in range(a+1,12)]
pools=[]
for fn in sys.argv[2:]:
    L=open(fn).read().split('\n'); n=int(L[0]); pools.append([list(map(int,L[1+i].split())) for i in range(n)])
c=Counter(); bad=[]
for ln in open(sys.argv[1]):
    left,right=ln.split('|')
    gens=[]
    for e in map(int,left.split()):
        z=[0]*12; z[E[e][0]]=1; z[E[e][1]]=-1; gens.append(z)
    for tok in right.split():
        p,i=map(int,tok.split(':')); gens.append(pools[p][i])
    Z=np.array([g[1:] for g in gens],dtype=np.int64)
    assert all(sum(g)==0 for g in gens)
    N=((Z@D.T)!=0).sum(0); assert (N<=tau).all()
    S=D[N==tau]
    r=nmod_mat(S.tolist(),PR).rank() if len(S) else 0
    c[r]+=1
    if r<11: bad.append(ln.strip())
print(sys.argv[1],'n',sum(c.values()),'rank dist',dict(c))
for b in bad[:20]: print('LOW RANK',b)
