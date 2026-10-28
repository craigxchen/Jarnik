import json, itertools, sys, math
import numpy as np
sys.path.insert(0,'.')
from rcore import EDGES, EIDX
perms=list(itertools.permutations(range(7)))
EP=np.array([[EIDX[tuple(sorted((s[i],s[j])))] for (i,j) in EDGES] for s in perms])
POW=np.array([7**k for k in range(21)],dtype=np.int64); Wm=POW[EP].T
AU='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/search_checks/'
for D in (4,5,6):
    au=json.load(open(AU+'reps_D%d.json'%D)); me=json.load(open('rreps_D%d.json'%D))
    def codes(reps):
        V=np.zeros((len(reps),21),dtype=np.int64)
        for r,A in enumerate(reps):
            for e in A: V[r,EIDX[tuple(e)]]+=1
        return set((V@Wm).min(axis=1).tolist()), len(reps)
    ca,na=codes(au); cm,nm=codes(me)
    # author's chunk coverage
    tot=0; seen={}
    for l in open(AU+'out/D%d.jsonl'%D):
        r=json.loads(l); seen.setdefault(r['orbit'],[]).append((r['c0'],r['c1'],r['n']))
    for o,iv in seen.items():
        iv.sort(); 
        assert iv[0][0]==0 and all(iv[k][1]==iv[k+1][0] for k in range(len(iv)-1)), (D,o,iv)
        tot+=sum(x[2] for x in iv)
        s=len(set(tuple(e) for e in au[o])); assert iv[-1][1]>=math.comb(21-s+D-1,D) or True
        last=iv[-1]; assert last[0]+last[2]==math.comb(21-s+D-1,D),(D,o,last)
    print(D,'author reps',na,'distinct orbits',len(ca),'equal to mine',ca==cm,'orbits covered',len(seen),'pairs in chunks',tot)
