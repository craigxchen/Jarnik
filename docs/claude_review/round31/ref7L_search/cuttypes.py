# enumerate 5-point layer types (squarefree: columns = cuts) with all pair distances >= h
import itertools, sys, pickle
from collections import Counter
PTS=range(5)
def canon_cut(S):
    S=frozenset(S)
    if 0 in S: S=frozenset(PTS)-S
    return S
CUTS=[]
for r in range(1,5):
    for S in itertools.combinations(range(1,5),r):
        CUTS.append(frozenset(S))
# include cuts as subsets not containing 0 (complement-canonical): sizes 1..4 of {1..4}
CUTS=[c for c in CUTS]
assert len(CUTS)==15
idx={c:i for i,c in enumerate(CUTS)}
PAIRS=list(itertools.combinations(PTS,2))
sep=[[1 if ((i in c)!=(j in c)) else 0 for (i,j) in PAIRS] for c in CUTS]
perms=list(itertools.permutations(PTS))
permmap=[[idx[canon_cut({p[x] for x in c})] for c in CUTS] for p in perms]
def canon(counts):
    best=None
    for pm in permmap:
        v=[0]*15
        for k,cnt in enumerate(counts):
            if cnt: v[pm[k]]+=cnt
        t=tuple(v)
        if best is None or t>best: best=t
    return best
def enum(D,h):
    res=set()
    maxsep=[max(sep[k][p] for k in range(15)) for p in range(10)]
    def rec(k,rem,counts,dist):
        if k==15:
            if rem==0 and all(d>=h for d in dist):
                res.add(canon(counts))
            return
        # prune: each column separates at most 1 per pair; need dist+rem>=h
        if any(d+rem<h for d in dist): return
        for c in range(rem,-1,-1):
            nd=[dist[p]+c*sep[k][p] for p in range(10)]
            counts.append(c)
            rec(k+1,rem-c,counts,nd)
            counts.pop()
    rec(0,D,[],[0]*10)
    return sorted(res)
if __name__=="__main__":
    out={}
    for D,h in [(7,4),(9,5),(10,6),(11,6),(12,7),(13,7)]:
        T=enum(D,h)
        out[(D,h)]=T
        print(D,h,len(T),flush=True)
    pickle.dump(out,open('sqfree_types.pkl','wb'))
