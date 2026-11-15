# Judge check (unrefereed). D=8 (or D given), c=0: all codes containing an orbit representative v0 of S_D x {+-1},
# pairwise balanced, pruned by saturated Lemma 1.1 (test box [-B,B]^D) and Gordan is implied for integer nonneg vectors.
# Outputs: max surviving code size, and canonical classes of surviving 6-point rank-5 codes.
import sys, itertools, numpy as np, flint, time, pickle
D=int(sys.argv[1]); c=int(sys.argv[2]); B=int(sys.argv[3]); MAXDEPTH=int(sys.argv[4])
K=D//2+c; L=D-1-K
pts=np.array(list(itertools.product([1,-1],repeat=D)),dtype=np.int64); n=len(pts)
dl=(pts[:,None,:]-pts[None,:,:])//2
p=(dl>0).sum(2); q=(dl<0).sum(2)
A=(p+q>=L+1)&(2*np.minimum(p,q)>=L)
adj=[set(np.nonzero(A[i])[0].tolist()) for i in range(n)]
T=[]
for a in itertools.product(range(-B,B+1),repeat=D):
    a=np.array(a)
    if not a.any(): continue
    ap=a[a>0].sum(); am=-a[a<0].sum()
    if ap+am<=L or 2*min(ap,am)<L: T.append(a)
T=np.array(T,dtype=np.int64).T
def killed(diffs):
    M=flint.fmpz_mat([[int(x) for x in d] for d in diffs])
    Y,nul=M.nullspace()
    if nul==0: return True
    Ym=np.array([[int(Y[r,cc]) for cc in range(nul)] for r in range(D)],dtype=np.int64)
    return bool((~(Ym.T@T).any(axis=0)).any())
perms6=np.array(list(itertools.permutations(range(6))))
def canon(rows):
    X=np.array(rows)  # 6 x D
    best=None
    for s in (1,-1):
        Y=(s*X>0).astype(np.int64)  # 6 x D bits
        Yp=Y[perms6]               # 720 x 6 x D
        codes=(Yp*(1<<np.arange(6))[None,:,None]).sum(1)  # 720 x D column codes
        codes.sort(axis=1)
        m=min(map(tuple,codes))
        if best is None or m<best: best=m
    return best
best=0; classes={}; t0=time.time(); nodes=0
def dfs(cl,cand,v0):
    global best,nodes
    nodes+=1
    if len(cl)>best: best=len(cl); print("best",best,time.time()-t0,flush=True)
    if len(cl)==6:
        diffs=[pts[x]-pts[v0] for x in cl[1:]]
        r=flint.fmpz_mat([[int(x) for x in d] for d in diffs]).rank()
        if r==5:
            k=canon([pts[x] for x in cl])
            if k not in classes: classes[k]=[list(map(int,pts[x])) for x in cl]
    if len(cl)>=MAXDEPTH: return
    last=cl[-1] if len(cl)>1 else -1
    for v in sorted(cand):
        if v<=last: continue
        new=cl+[v]
        if len(new)>=3 and killed([pts[x]-pts[v0] for x in new[1:]]): continue
        dfs(new,cand&adj[v],v0)
for w in range(0,D//2+1):
    v0=int(np.nonzero(((pts>0).sum(1)==w)&(np.all(pts[:, :w]>0,axis=1)))[0][0])
    dfs([v0],adj[v0],v0)
    print("rep weight",w,"done; best",best,"classes",len(classes),"nodes",nodes,"t",time.time()-t0,flush=True)
pickle.dump(classes,open(f"classes_D{D}_c{c}.pkl","wb"))
print("FINAL best",best,"rank5 6-pt classes",len(classes))
for k,v in classes.items(): print(v)
