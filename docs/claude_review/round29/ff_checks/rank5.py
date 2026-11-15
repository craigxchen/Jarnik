import itertools, numpy as np
from zerosum import solve, count_zero_sum
# extended Hamming [8,4,4]
G=np.array([[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]])
H8=sorted({tuple(np.array(c)@G%2) for c in itertools.product([0,1],repeat=4)})
allv=np.array(list(itertools.product([0,1],repeat=8)))
H8nz=[np.array(h) for h in H8 if any(h)]
seen=set(); cands=[]
for comb in itertools.combinations(range(len(H8nz)),5):
    B=np.array([H8nz[i] for i in comb],float)
    if np.linalg.matrix_rank(B)<5: continue
    # 0/1 vectors in span
    P=B.T@np.linalg.pinv(B.T)
    inU=np.all(np.abs(allv@P.T-allv)<1e-9,axis=1)
    code=frozenset(map(tuple,allv[inU]))
    if code in seen: continue
    seen.add(code)
    w=[sum(c) for c in code if any(c)]
    if min(w)<4: continue
    if len(code)>=9: cands.append((len(code),code,B))
print("distinct rank-5 spans with >=9 words, all weights>=4:",len(cands), sorted(set(c[0] for c in cands)))
# test realizability (up to coordinate permutation equivalence we just test all; small)
import collections
res=collections.Counter()
for n,code,B in cands[:40]:
    Fs=[list(np.nonzero(b)[0]) for b in B.astype(int)]
    sols=solve(Fs,3,8,tries=40)
    nd=[d for d,r in sols if d>1e-3]
    res[(n,len(nd)>0)]+=1
    if nd: print("REALIZED",n,code)
print(res)
