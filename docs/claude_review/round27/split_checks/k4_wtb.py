# T2 control at k=4: near-tight split profile polytope Pi_nt^(4) and the k=4 analogue of WTB: 2 reg(S) <= cost_x(S).
import itertools, numpy as np, random, sys, time
from flint import nmod_mat
from scipy.optimize import linprog
n=4
Bs=[frozenset(c) for r in range(1,n) for c in itertools.combinations(range(n),r)]
idx={B:i for i,B in enumerate(Bs)}; comp=lambda B: frozenset(range(n))-B
A_eq=[];b_eq=[];seen=set()
for B in Bs:
    C=comp(B)
    if C in seen: continue
    seen.add(B); row=[0]*14; row[idx[B]]=1; row[idx[C]]=1; A_eq.append(row); b_eq.append(1)
for a in range(1,n):
    row=[0]*14
    for B in Bs:
        if 0 in B: row[idx[B]]+=1
        if a in B: row[idx[B]]-=1
    A_eq.append(row); b_eq.append(0)
A_eq=np.array(A_eq,float); b_eq=np.array(b_eq,float)
import sympy as sp
print("dim Pi_nt^(4) =",14-sp.Matrix(A_eq).rank())
Bmat=np.array([[1 if a in B else 0 for a in range(n)] for B in Bs])
def full(m): return tuple(list(m)+[-sum(m)])
def pis(S): S4=np.array([full(m) for m in S]); return (-(S4@Bmat.T)).max(axis=0)
def monos(d): return [a for t in range(d+1) for a in itertools.product(range(t+1),repeat=3) if sum(a)==t]
def reg(S,p=2**31-1):
    S=list(S)
    for d in range(0,60):
        ms=monos(d)
        if len(ms)<len(S): continue
        rows=[[pow(m[0]%p,a[0],p)*pow(m[1]%p,a[1],p)%p*pow(m[2]%p,a[2],p)%p for a in ms] for m in S]
        if nmod_mat(rows,p).rank()==len(S): return d
def mincost(S):
    r=linprog(pis(S).astype(float),A_eq=A_eq,b_eq=b_eq,bounds=[(0,1)]*14,method='highs'); return r.fun
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 0)
best=0;t0=time.time();cnt=0
while time.time()-t0<120:
    # random lattice points in a random polytope {<m,1_B> >= -p_B}
    p=[random.randint(0,3) for _ in Bs]
    pts=[m for m in itertools.product(range(-6,7),repeat=3) if all(-sum(full(m)[a] for a in B)<=p[idx[B]] for B in Bs)]
    if len(pts)<2 or len(pts)>600: continue
    if random.random()<0.5:
        pts=random.sample(pts,max(2,len(pts)*2//3))
    r=reg(pts); mc=mincost(pts); cnt+=1
    if mc>0 and 2*r/mc>best-1e-12:
        best=max(best,2*r/mc); print(cnt,"|S|",len(pts),"reg",r,"mincost",round(mc,4),"ratio 2reg/mincost",round(2*r/mc,4),flush=True)
print("done",cnt,"best",best)
