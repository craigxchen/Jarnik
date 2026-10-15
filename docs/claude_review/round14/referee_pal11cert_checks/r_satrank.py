# saturated-rank test for a list of 7-multisets of generators (12-vectors)
import numpy as np, pickle, sys
from flint import nmod_mat
d=pickle.load(open('r_box.pkl','rb')); D=d['D'].astype(np.int64); tau=d['tau']
PR=2147483647
def satrank(gens):
    Z=np.array([g[1:] for g in gens],dtype=np.int64)
    N=((Z@D.T)!=0).sum(0)
    assert (N<=tau).all()
    S=D[N==tau]
    if len(S)==0: return 0,0
    r=nmod_mat(S.tolist(),PR).rank()
    return r,len(S)
if __name__=='__main__':
    E=[(a,b) for a in range(12) for b in range(a+1,12)]
    def pairvec(e):
        z=[0]*12; z[E[e][0]]=1; z[E[e][1]]=-1; return z
    from collections import Counter
    c=Counter()
    for ln in open(sys.argv[1]):
        es=list(map(int,ln.split()))
        r,ns=satrank([pairvec(e) for e in es]); c[r]+=1
    print(sys.argv[1],'saturated rank distribution',dict(c))
