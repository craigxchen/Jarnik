import pickle, itertools, numpy as np
from flint import nmod_mat
d=pickle.load(open('r_box.pkl','rb')); A=d['A'].astype(int)
def passes(S):
    B=np.unique(A[:,list(S)],axis=0); s=len(S)
    mons=[a for a in itertools.product((0,1,2),repeat=s) if sum(a)<=7]
    K=[[int(np.prod([b[i]**a[i] for i in range(s)])) for a in mons] for b in B]
    return nmod_mat(K,2147483647).rank()==len(B)
P4=[S for S in itertools.combinations(range(11),4) if passes(S)]
P5=[S for S in itertools.combinations(range(11),5) if passes(S)]
maxi4=[S for S in P4 if not any(set(S)<=set(T) for T in P5)]
print(len(P4),len(P5),'maximal 4-sets',len(maxi4))
