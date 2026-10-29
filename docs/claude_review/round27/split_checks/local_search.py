# Local search for S violating Conjecture WTB: maximise 4*reg(S) - min_{x in Pi_nt} cost_x(S).
# Start from root zonotopes (Theorem-B tight) and perturb by deleting/adding lattice points.
import sys, time, random; sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks')
from reglib import *
from scipy.optimize import linprog
A_eq=[];b_eq=[];seen=set()
for B in Bs:
    C=comp(B)
    if C in seen: continue
    seen.add(B); row=[0]*30; row[idx[B]]=1; row[idx[C]]=1; A_eq.append(row); b_eq.append(1)
for a in range(1,5):
    row=[0]*30
    for B in Bs:
        if 0 in B: row[idx[B]]+=1
        if a in B: row[idx[B]]-=1
    A_eq.append(row); b_eq.append(0)
A_eq=np.array(A_eq,float); b_eq=np.array(b_eq,float)
def mincost(S):
    pi=pis(S); r=linprog(pi.astype(float),A_eq=A_eq,b_eq=b_eq,bounds=[(0,1)]*30,method='highs'); return r.fun
roots=[]
for a in range(5):
    for b in range(5):
        if a<b:
            v=[0]*5; v[a]=1; v[b]=-1; roots.append(tuple(v[:4]))
def zonotope(gens):
    pts={(0,0,0,0)}
    for g,k in gens:
        new=set()
        for p in pts:
            for j in range(k+1): new.add(tuple(p[i]+j*g[i] for i in range(4)))
        pts=new
    return pts
random.seed(int(sys.argv[1])); t0=time.time(); T=float(sys.argv[2])
best=-1e9
while time.time()-t0<T:
    gens=[(random.choice(roots),random.randint(1,2)) for _ in range(random.randint(2,5))]
    S=set(zonotope(gens))
    if len(S)>250: continue
    def score(S):
        r=reg(list(S)); return 4*r-mincost(list(S)), r
    cur,rr=score(S)
    for step in range(60):
        S2=set(S)
        if random.random()<0.6 and len(S2)>2:
            S2.discard(random.choice(list(S2)))
        else:
            p=random.choice(list(S2)); g=random.choice(roots); s=random.choice([1,-1])
            S2.add(tuple(p[i]+s*g[i] for i in range(4)))
        if len(S2)>250: continue
        sc,r2=score(S2)
        if sc>=cur-1e-9:
            S,cur,rr=S2,sc,r2
    if cur>best-1e-9:
        best=max(best,cur); print(f"best 4reg-mincost={cur:.4f} reg={rr} |S|={len(S)}",flush=True)
print("done best",best)
