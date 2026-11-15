# Enumerate equidistant (pairwise orthogonal) +-1 codes containing 1, length D, size M, up to coordinate permutation
# (column-type multiset). Output type-count vectors.
import itertools, sys, json
from collections import Counter
def classes(D,M):
    # a code with s1=1 is given by a multiset of column types t in {+-1}^{M-1} (column k = (s_2k..s_Mk)); orthogonality:
    # for i: sum_t n_t t_i = 0 ; for i<j: sum_t n_t t_i t_j = 0; sum n_t = D. Column-sign flips change s1 so we keep s1=1
    # but we may also flip a column's sign together with... (flip changes s1_k) -> not allowed while s1=1 fixed; fine.
    T=list(itertools.product([1,-1],repeat=M-1))
    res=[]
    def rec(idx,remaining,cnt):
        if idx==len(T):
            if remaining: return
            ok=all(sum(n*t[i] for t,n in cnt.items())==0 for i in range(M-1)) and \
               all(sum(n*t[i]*t[j] for t,n in cnt.items())==0 for i in range(M-1) for j in range(i+1,M-1))
            if ok: res.append(dict(cnt))
            return
        for n in range(remaining+1):
            if n: cnt[T[idx]]=n
            rec(idx+1,remaining-n,cnt)
            if n: del cnt[T[idx]]
    rec(0,D,{})
    return res
if __name__=='__main__':
    D=int(sys.argv[1]);M=int(sys.argv[2])
    r=classes(D,M)
    # reduce by permutations of points 2..M and by distinctness (columns may repeat -> repeated types allowed; roots distinct anyway)
    canon=set()
    out=[]
    for c in r:
        best=None
        for perm in itertools.permutations(range(M-1)):
            key=tuple(sorted((tuple(t[p] for p in perm),n) for t,n in c.items()))
            if best is None or key<best: best=key
        if best not in canon:
            canon.add(best); out.append(best)
    print(len(r),'raw',len(out),'classes')
    json.dump(out,open('equi_%d_%d.json'%(D,M),'w'))
