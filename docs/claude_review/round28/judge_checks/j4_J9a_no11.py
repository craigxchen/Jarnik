# Judge check J4 (unrefereed): independent proof that syl16|J9a (L=3) contains no lattice zonotope of order 11
# (any primitive generators in {-2..2}^8). Own LP dual (rationalised, checked exactly over all candidates),
# then exhaustive DFS over the generators allowed by the dual.
import itertools, math, sys, time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
def dot(u,x): return bin(u&x).count('1')&1
cuts=set()
J=list(range(9)); Js=set(J)
for u in range(1,16):
    S=frozenset(x for x in range(16) if dot(u,x))&Js
    if S and S!=Js: cuts.add(min(S,frozenset(Js-S),key=lambda s:(len(s),sorted(s))))
cuts=sorted(cuts,key=sorted); L=3; d=8
A=[k for k in itertools.product((-1,0,1),repeat=d) if all(abs(sum((list(k)+[-sum(k)])[x] for x in S))<=L for S in cuts)]
print('|A|',len(A),'cuts',len(cuts)); Aset=set(A)
# functionals: coordinates, base coordinate sigma, slabs (as vectors on k')
lams=[]
for x in range(d): v=[0]*d; v[x]=1; lams.append(v)
lams.append([1]*d)
for S in cuts:
    v=[(1 if x in S else 0)-(1 if 8 in S else 0) for x in range(d)]; lams.append(v)
Am=np.array(A); Lm=np.array(lams)
vals=Am@Lm.T; W=vals.max(0)-vals.min(0)
print('widths',W.tolist())
cand=[]
for v in itertools.product(range(-2,3),repeat=d):
    nz=[c for c in v if c]
    if not nz or nz[0]<0: continue
    if math.gcd(*[abs(c) for c in v])!=1: continue
    cand.append(v)
C=np.array(cand); absv=np.abs(C@Lm.T)
ok=(absv<=W).all(1); C=C[ok]; absv=absv[ok]
print('candidates fitting widths',len(C))
res=linprog(W.astype(float),A_ub=-absv.astype(float),b_ub=-np.ones(len(C)),bounds=[(0,None)]*len(W),method='highs')
print('LP float',res.fun)
y=[F(float(t)).limit_denominator(1000) for t in res.x]
den=math.lcm(*[t.denominator for t in y]); yi=[int(t*den) for t in y]   # integers
cv=absv@np.array(yi,dtype=np.int64)              # exact integer c_v * den
m=int(cv.min()); B=F(int(np.dot(W,yi)),m)        # rescale so min c_v = 1
print('exact dual bound B =',B,'(min c_v scaled to 1)')
t=11
allowed=[tuple(C[i]) for i in range(len(C)) if F(int(cv[i]),m)-1<=B-t]
print('generators allowed for order',t,':',len(allowed))
# DFS on boolean grid
idx=lambda k: tuple(c+1 for c in k)
grid=np.zeros((3,)*d,bool)
for k in A: grid[idx(k)]=True
def shift_and(K,v):
    # K ∩ (K - v): k in new iff k in K and k+v in K
    out=np.zeros_like(K); src=[];dst=[]
    for c in v:
        if c>=0: dst.append(slice(0,3-c)); src.append(slice(c,3))
        else: dst.append(slice(-c,3)); src.append(slice(0,3+c))
    if any(s.start>=s.stop for s in src): return out
    out[tuple(dst)]=K[tuple(dst)]&K[tuple(src)]; return out
nodes=0; found=None; t0=time.time()
def dfs(K,start,depth,path):
    global nodes,found
    nodes+=1
    if depth==t: found=list(path); return True
    for i in range(start,len(allowed)):
        K2=shift_and(K,allowed[i])
        if K2.any():
            path.append(allowed[i])
            if dfs(K2,i,depth+1,path): return True
            path.pop()
    return False
r=dfs(grid,0,0,[])
print('order-11 zonotope exists:',r,'nodes',nodes,'secs',round(time.time()-t0))
# sanity: same search finds order 10
t=10; nodes=0
allowed10=[tuple(C[i]) for i in range(len(C)) if F(int(cv[i]),m)-1<=B-t]
allowed=allowed10
r10=dfs(grid,0,0,[]); print('order-10 found:',r10,'gens',found,'nodes',nodes)
