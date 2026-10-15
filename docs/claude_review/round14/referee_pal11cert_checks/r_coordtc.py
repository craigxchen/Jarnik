# Referee: which coordinate sets S have rank_Q K_{<=7}(pi_S A) = |pi_S A| (=> span(e_S) in TC(A,8) via fatpoint Prop 2.2(b), Lemma 2.1).
import pickle, itertools, numpy as np
from flint import fmpz_mat, nmod_mat
d=pickle.load(open('r_box.pkl','rb')); A=d['A'].astype(int); twin=d['twin']; core=d['core']
def monos(s,deg):
    return [a for a in itertools.product((0,1,2),repeat=s) if sum(a)<=deg]
res={}
for size in (3,4,5,6):
    passed=[]
    for S in itertools.combinations(range(11),size):
        B=np.unique(A[:,list(S)],axis=0)
        mons=monos(size,7)
        K=np.array([[int(np.prod([b[i]**a[i] for i in range(size)])) for a in mons] for b in B])
        r=nmod_mat(K.tolist(),2147483647).rank()
        if r==len(B): passed.append(S)
        elif size==6:
            rq=fmpz_mat(K.tolist()).rank()
            assert rq<len(B)
    res[size]=passed
    print('size',size,'passing sets',len(passed))
tw=set(tuple(sorted(x-1 for x in P)) for P in twin)
print('5-sets passing == twin patterns:', set(res[5])==tw)
