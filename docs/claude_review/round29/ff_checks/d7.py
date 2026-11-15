# D=7, c=0 (K=3, L=3): try to realise 5- and 4-point clusters (zero-sum families) on the moment surface.
import itertools, numpy as np
from zerosum import solve, count_zero_sum
D=7;L=3
words=[w for w in itertools.product([0,1],repeat=D) if sum(w)>=4]
def dist(a,b): return sum(x!=y for x,y in zip(a,b))
allv=np.array(list(itertools.product([0,1],repeat=D)),float)
done=set(); found={}
for k in (4,3):
    for comb in itertools.combinations(words,k):
        if any(dist(a,b)<4 for a,b in itertools.combinations(comb,2)): continue
        if comb[0]!=(1,1,1,1,0,0,0): continue   # symmetry normalisation (S_7 acts transitively on weight-4 words); weight>=5 first words skipped
        B=np.array(comb,float)
        if np.linalg.matrix_rank(B)<k: continue
        P=B.T@np.linalg.pinv(B.T); inU=np.all(np.abs(allv@P.T-allv)<1e-9,axis=1)
        code=allv[inU]
        if any(dist(a,b)<4 for a,b in itertools.combinations(map(tuple,code),2)): continue
        key=tuple(sorted(map(lambda r:tuple(int(x) for x in r),code)))
        if key in done: continue
        done.add(key)
        Fs=[list(np.nonzero(b)[0]) for b in B.astype(int)]
        sols=solve(Fs,L,D,tries=30)
        nd=[(d,r) for d,r in sols if d>1e-3]
        if nd:
            d,r=max(nd,key=lambda s:s[0]); z=count_zero_sum(r,L)
            print("REALIZED: code size",len(code),"rank",k,"-> zero-sum subsets found",len(z),"dmin",round(d,4),flush=True)
            found[len(code)]=found.get(len(code),0)+1
    print("k",k,"done; realized sizes:",found,flush=True)
