import itertools, numpy as np, sys
from zerosum import solve, count_zero_sum
D=int(sys.argv[1]);L=int(sys.argv[2]); ks=[int(x) for x in sys.argv[3].split(',')]; tries=int(sys.argv[4])
words=[w for w in itertools.product([0,1],repeat=D) if sum(w)>=L+1]
def dist(a,b): return sum(x!=y for x,y in zip(a,b))
allv=np.array(list(itertools.product([0,1],repeat=D)),float)
perms=list(itertools.permutations(range(D)))
def canon(code):
    return min(tuple(sorted(tuple(w[p] for p in pm) for w in code)) for pm in perms)
done=set()
for k in ks:
    found={}
    for comb in itertools.combinations(words,k):
        if any(dist(a,b)<L+1 for a,b in itertools.combinations(comb,2)): continue
        B=np.array(comb,float)
        if np.linalg.matrix_rank(B)<k: continue
        P=B.T@np.linalg.pinv(B.T); inU=np.all(np.abs(allv@P.T-allv)<1e-9,axis=1)
        code=[tuple(int(x) for x in r) for r in allv[inU]]
        if any(dist(a,b)<L+1 for a,b in itertools.combinations(code,2)): continue
        key=canon(code)
        if key in done: continue
        done.add(key)
        Fs=[list(np.nonzero(b)[0]) for b in B.astype(int)]
        sols=solve(Fs,L,D,tries=tries)
        nd=[(d,r) for d,r in sols if d>1e-3]
        tag=len(code)
        if nd:
            d,r=max(nd,key=lambda s:s[0]); z=count_zero_sum(r,L)
            found[tag]=found.get(tag,0)+1
            print("REALIZED rank",k,"code size",tag,"actual cluster",len(z),"dmin",round(d,4),key,flush=True)
        else: print("not realized rank",k,"size",tag,key,flush=True)
    print("k",k,"done",found,flush=True)
