# Independent enumeration of pair-admissible 6-row sign codes at D=8,c=0 up to S_8 x global flip.
import itertools, numpy as np, pickle, sys
D=8
rows=list(range(256))
def bits(r): return [(r>>j)&1 for j in range(D)]
def adm(r,s):
    p=bin(r&~s&255).count('1'); n=bin(s&~r&255).count('1')
    return p>=2 and n>=2
adj={r:set(s for s in rows if s!=r and adm(r,s)) for r in rows}
# permutation tables
perms=list(itertools.permutations(range(D)))
P=np.array(perms,dtype=np.int64)          # 40320 x 8
B=np.array([bits(r) for r in rows],dtype=np.int64)  # 256x8
# image of row r under perm p: bit p[j] of image = bit j of r   (column j -> position p[j])
T=np.zeros((len(perms),256),dtype=np.int64)
for j in range(D):
    T+= (B[:,j][None,:]) << P[:,j][:,None]
T=T.astype(np.uint16)
def canon(S):
    S=list(S)
    best=None
    for fl in (0,255):
        img=T[:,[r^fl for r in S]]           # 40320 x k
        img=np.sort(img,axis=1)
        # lexicographic min
        packed=np.zeros(img.shape[0],dtype=np.int64)
        for c in range(img.shape[1]): packed=(packed<<8)|img[:,c].astype(np.int64)
        i=int(np.argmin(packed)); key=tuple(int(v) for v in img[i])
        if best is None or key<best: best=key
    return best
# enumerate all cliques of size k by extension of canonical reps (orbit BFS), complete since every
# k-set restricts to a (k-1)-set (any row removal), which is in some orbit; extending the rep by all rows
# and canonising gives every k-orbit.
level={(): None}
reps=[()]
counts=[]
for k in range(1,7):
    new={}
    for R in reps:
        cand=set(rows) if not R else set.intersection(*[adj[r] for r in R])
        for s in cand:
            if s in R: continue
            c=canon(R+(s,))
            if c not in new: new[c]=1
    reps=sorted(new)
    counts.append(len(reps)); print(k,len(reps),flush=True)
    pickle.dump(reps,open(f'reps{k}.pkl','wb'))
