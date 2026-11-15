# Judge check (unrefereed): exhaustive DFS over codes in {+-1}^D at level c, pairwise balanced,
# pruned by the saturated-lattice lemma (Lemma 1.1 search = Lemma 1.2 prove):
#   every nonzero integer a in span_Q(differences) has ||a||_1 >= L+1 and min(a+,a-) >= L/2.
# Kills are exact (flint integer nullspace, int64 tests). Test vectors: all a with ||a||_1<=L, and all
# a in [-B,B]^D with min(a+,a-) < L/2 (sound subset of the full condition).
import sys, itertools, numpy as np, flint, time
D=int(sys.argv[1]); c=int(sys.argv[2]); B=int(sys.argv[3]) if len(sys.argv)>3 else 2
K=D//2+c; L=D-1-K
pts=np.array(list(itertools.product([1,-1],repeat=D)),dtype=np.int64)
n=len(pts)
def bal(i,j):
    d=(pts[i]-pts[j])//2; p=(d>0).sum(); q=(d<0).sum()
    return p+q>=L+1 and 2*min(p,q)>=L
adj=[set() for _ in range(n)]
for i in range(n):
    for j in range(i+1,n):
        if bal(i,j): adj[i].add(j); adj[j].add(i)
# test vectors
T=[]
for a in itertools.product(range(-B,B+1),repeat=D):
    a=np.array(a); 
    if not a.any(): continue
    ap=a[a>0].sum(); am=-a[a<0].sum()
    if ap+am<=L or 2*min(ap,am)<L: T.append(a)
# also weight<=L vectors with entries beyond B (only if L>B)
if L>B:
    for a in itertools.product(range(-L,L+1),repeat=D):
        a=np.array(a)
        if a.any() and np.abs(a).sum()<=L and np.abs(a).max()>B: T.append(a)
T=np.array(T,dtype=np.int64).T  # D x nT
print("D",D,"c",c,"L",L,"edges",sum(len(s) for s in adj)//2,"tests",T.shape[1],flush=True)
def killed(diffs):
    M=flint.fmpz_mat([[int(x) for x in d] for d in diffs])
    # complement: vectors y with M y = 0 ; a in rowspace(M) iff Y^T a = 0 for basis Y of nullspace(M)
    Y,nul=M.nullspace()
    if nul==0: return True  # span is everything -> contains e_1 (weight 1)
    Ym=np.array([[int(Y[r,cc]) for cc in range(nul)] for r in range(D)],dtype=object).astype(np.int64)
    R=Ym.T@T
    return bool((~R.any(axis=0)).any())
def rank(diffs):
    return flint.fmpz_mat([[int(x) for x in d] for d in diffs]).rank()
best=0; stats={}; cnt=0; t0=time.time()
seen_rank5=[]
def dfs(cl,cand):
    global best,cnt
    cnt+=1
    m=len(cl)
    if m>=2:
        r=rank([pts[x]-pts[cl[0]] for x in cl[1:]])
        stats[(m,r)]=stats.get((m,r),0)+1
        if m>best: best=m; print("new best",m,"rank",r,[list(pts[x]) for x in cl],flush=True)
    for v in sorted(cand):
        if v<=cl[-1]: continue
        new=cl+[v]
        diffs=[pts[x]-pts[new[0]] for x in new[1:]]
        if len(new)>=3 and killed(diffs): continue
        dfs(new,cand&adj[v])
# symmetry: S_D acts; only need clique whose min-index vertex ... use orbit reps for first vertex is not valid
# for increasing-order enumeration, so enumerate all.
for v0 in range(n):
    dfs([v0],adj[v0])
print("best",best,"nodes",cnt,"time",time.time()-t0)
print(sorted(stats.items()))
