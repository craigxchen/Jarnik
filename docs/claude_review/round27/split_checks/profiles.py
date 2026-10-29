# Near-tight split profile polytope Pi_nt and binomial section costs c_x(m) = sum_B x_B (m(B))^+.
# Pi_nt = {x in [0,1]^30 : x_B + x_{B^c} = 1, P_ab(x) := sum_{B ni a, not b} x_B = 4 for all ordered a != b}.
import itertools, numpy as np, sympy as sp
from fractions import Fraction
from scipy.optimize import linprog
n=5
Bs=[frozenset(c) for r in range(1,n) for c in itertools.combinations(range(n),r)]
idx={B:i for i,B in enumerate(Bs)}
comp=lambda B: frozenset(range(n))-B
A_eq=[];b_eq=[]
for B in Bs:
    if min(B)==0 and True:
        pass
done=set()
for B in Bs:
    C=comp(B)
    if (C,B) in done: continue
    done.add((B,C))
    row=[0]*30; row[idx[B]]=1; row[idx[C]]=1; A_eq.append(row); b_eq.append(1)
for a in range(n):
    for b in range(n):
        if a==b: continue
        row=[0]*30
        for B in Bs:
            if a in B and b not in B: row[idx[B]]=1
        A_eq.append(row); b_eq.append(4)
M=sp.Matrix(A_eq)
print("rank of equality system:", M.rank(), " => dim Pi_nt =", 30-M.rank())
# check x=1/2 feasible
x0=np.full(30,0.5); print("x=1/2 residual", np.abs(np.array(A_eq)@x0-np.array(b_eq)).max())
# W profile (labels 0..4) from round 23: F list
F=[{0,3,4},{0,1,3,4},{3,4},{1,3,4},{0,2,4},{0,1,2,4},{2,4},{1,2,4},{0,2,3},{0,1,2,3},{2,3},{1,2,3},{1},{0},{0,1}]
xW=np.zeros(30)
for B in F: xW[idx[frozenset(B)]]=1
print("x_W residual", np.abs(np.array(A_eq)@xW-np.array(b_eq)).max())
# vertex enumeration by random LPs
rng=np.random.default_rng(0)
verts={}
for it in range(4000):
    c=rng.normal(size=30)
    r=linprog(c,A_eq=np.array(A_eq,float),b_eq=np.array(b_eq,float),bounds=[(0,1)]*30,method='highs')
    if r.status==0:
        v=tuple(np.round(r.x,6))
        verts[v]=verts.get(v,0)+1
print("distinct vertices found:", len(verts))
V=np.array(list(verts.keys()))
vals=sorted(set(np.round(V.flatten(),6)))
print("coordinate values at vertices:", vals[:20])
nint=sum(1 for v in V if all(abs(t-round(t))<1e-9 for t in v))
print("0/1 vertices:", nint)
# binomial costs
ms=[m for m in itertools.product(range(-3,4),repeat=5) if sum(m)==0 and any(m)]
def mB(m,B): return sum(m[a] for a in B)
Cmat=np.array([[max(0,mB(m,B)) for B in Bs] for m in ms],float)
best=[]
for v in V:
    cost=Cmat@v
    i=np.argmin(cost); best.append((cost[i],ms[i]))
mins=sorted(set(round(b[0],6) for b in best))
print("min binomial cost over vertices: distinct values", mins[:10])
costW=Cmat@xW; print("min cost at x_W:", costW.min(), "argmin examples", [ms[i] for i in np.where(costW<=costW.min()+1e-9)[0][:6]])
cost05=Cmat@x0; print("min cost at 1/2:", cost05.min())
np.save('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/pi_nt_vertices.npy',V)
