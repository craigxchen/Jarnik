import numpy as np, itertools, math, time
from flint import nmod_mat
from reng import REng, PR
from rlib import J_exact, monos
d=5; S=[{0,1,4},{0,2,3,4},{1,2,3,4}]; Ls=[1,1,1]; DT=int(__import__('sys').argv[1])
COLS=np.array([[1 if x in s else 0 for x in range(d)] for s in S],dtype=np.int64)
K=np.array(list(itertools.product([-1,0,1],repeat=d)),dtype=np.int64)
A=K[np.all(np.abs(K@COLS.T)<=np.array(Ls),axis=1)]; print("|A|",len(A))
E=REng(d,DT)
JA=J_exact([tuple(int(x) for x in p) for p in A],DT,cap=2)
print("dim J(A)_D exact",[len(JA[D]) for D in range(DT+1)])
class Acc:
    def __init__(s): s.B=None; s.r=0
    def add(s,rows):
        if not rows: return 0
        M=np.vstack(([s.B] if s.B is not None else [])+[np.array(rows)])
        nm=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR); R,rk=nm.rref()
        k=rk-s.r
        if k>0: s.B=np.array([[int(R[i,j]) for j in range(M.shape[1])] for i in range(rk)],dtype=np.int64); s.r=rk
        return k
def fam_rows(U):
    U=[np.array(u,dtype=np.int64) for u in U]
    img=np.unique(np.stack([A@u for u in U],axis=1),axis=0)
    sizes=[len(np.unique(img[:,i])) for i in range(len(U))]
    if np.prod(sizes)==len(img): return None
    J=J_exact([tuple(int(x) for x in p) for p in img],DT)
    rows=[]
    for D in range(DT+1):
        for vec in J[D]:
            g=E.pullback(vec,monos(len(U),D),U)
            if g.any(): rows+=E.ideal_rows(g,D,DT)
    return rows
dirs=[v for v in itertools.product(range(-2,3),repeat=d) if any(v) and math.gcd(*[abs(x) for x in v])==1 and next(x for x in v if x)>0]
acc=Acc()
for v in dirs:
    s=len(np.unique(A@np.array(v)))
    if s<=DT: acc.add(E.ideal_rows(E.pullback([1],[(s,)],[list(v)]),s,DT))
print("D",DT,"J_w(R=2)",acc.r,"J(A)",len(JA[DT]))
cols=[list(c) for c in COLS]
for F in itertools.combinations(range(3),3):
    r=fam_rows([cols[i] for i in F]); 
    if r: acc.add(r)
print("+ 3-col projection",acc.r)
d1=[v for v in dirs if max(abs(x) for x in v)<=1]
t0=time.time()
for U in itertools.combinations(d1,2):
    r=fam_rows(list(U))
    if r:
        k=acc.add(r)
        if k: print(" pair",U,"new",k,"total",acc.r)
print("after all pairs of {-1,0,1} dirs:",acc.r,"J(A)",len(JA[DT]),time.time()-t0)
