# Exact max cluster for 'free columns' model: columns arbitrary in R^L with positive first coordinate, distinct.
# For L<=2 this IS the true function-field problem (moment columns (y,2xy) are arbitrary with y>0).
# M = max |closure(U)|, U=span(S-s1), subject to: dim U >= D-L, Gordan (U has no nonzero x>=0), separation (e_k-e_l not in U).
import itertools, numpy as np, sys
from scipy.optimize import linprog
D=int(sys.argv[1]); L=int(sys.argv[2])
V=np.array(list(itertools.product([1,-1],repeat=D)),float)
def basis(vecs):
    if not vecs: return np.zeros((0,D))
    A=np.array(vecs); u,s,vt=np.linalg.svd(A); r=int((s>1e-9).sum()); return vt[:r]
def closure(s1,B):
    if B.shape[0]==0: return [tuple(s1)]
    P=B.T@B
    X=V-s1
    ok=np.all(np.abs(X@P-X)<1e-8,axis=1)
    return [tuple(v) for v in V[ok]]
def gordan(B):
    if B.shape[0]==0: return True
    # max sum x, x = B^T c, 0<=x<=1
    r=B.shape[0]
    res=linprog(-B.sum(axis=0) if False else -np.ones(D)@B.T, A_ub=np.vstack([B.T,-B.T]), b_ub=np.concatenate([np.ones(D),np.zeros(D)]), bounds=[(None,None)]*r)
    return res.status==0 and -res.fun<1e-7
def separation(B):
    if B.shape[0]==0: return True
    P=B.T@B
    for k,l in itertools.combinations(range(D),2):
        e=np.zeros(D); e[k]=1; e[l]=-1
        if np.linalg.norm(e@P-e)<1e-8: return False
    return True
best=1; seen=set(); bestS=None
# by symmetry (negation of all) s1 ranges over V; BFS over closures
frontier=[]
for s1 in V[:1] if False else V:
    frontier.append((s1,[]))
cnt=0
while frontier:
    s1,diffs=frontier.pop()
    B=basis(diffs)
    S=closure(s1,B)
    key=frozenset(S)
    if key in seen: continue
    seen.add(key)
    if not (gordan(B) and (L==1 or separation(B))): continue
    if len(S)>best: best=len(S); bestS=S; print("best",best,"dimU",B.shape[0],flush=True)
    if B.shape[0]>=D-1: continue
    Sset=set(S)
    for v in V:
        if tuple(v) in Sset: continue
        d=v-s1
        if np.sum(np.abs(d)>0)<L+1: continue
        frontier.append((s1,diffs+[d]))
    cnt+=1
print("D",D,"L",L,"max M",best, "closures examined",len(seen))
