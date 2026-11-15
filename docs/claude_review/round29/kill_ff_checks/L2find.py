# random closure growth for L=2 to find a large admissible closed set; then exact realisation with sympy
import itertools, numpy as np, sys, random, sympy as sp
D=int(sys.argv[1]); target=int(sys.argv[2]); K=int(sys.argv[3])
allv=np.array(list(itertools.product([0,1],repeat=D)),dtype=float)
bad=[np.eye(D)[k] for k in range(D)]
for k,l in itertools.combinations(range(D),2):
    bad.append(np.eye(D)[k]-np.eye(D)[l]); bad.append(np.eye(D)[k]+np.eye(D)[l])
bad=np.array(bad)
def closure(idx):
    B=allv[list(idx)]; u,s,vt=np.linalg.svd(B); rk=int((s>1e-9).sum()); Nn=vt[rk:]
    if Nn.shape[0]==0: return None
    if (np.abs(bad@Nn.T).max(axis=1)<1e-9).any(): return None
    return frozenset(np.nonzero(np.abs(allv@Nn.T).max(axis=1)<1e-9)[0].tolist())
random.seed(3); best=None
for trial in range(4000):
    S=frozenset([0])
    while True:
        cands=list(range(1,2**D)); random.shuffle(cands); grown=False
        for v in cands:
            if v in S: continue
            T=closure(S|{v})
            if T is not None: S=T; grown=True; break
        if not grown: break
    if best is None or len(S)>len(best): best=S; print("trial",trial,len(S),flush=True)
    if len(best)>=target: break
S=sorted(best); print([format(i,'0%db'%D) for i in S])
# exact realisation: W = 2 random integer combos of a rational basis of span(S)^perp
M=sp.Matrix([[int(c) for c in format(i,'0%db'%D)] for i in S])
perp=M.nullspace()
random.seed(7)
t=sp.symbols('t')
for attempt in range(20):
    rows=[sum((random.randint(-9,9)*v for v in perp),sp.zeros(D,1)) for _ in range(2)]
    W=sp.Matrix.hstack(*rows).T
    cols=[tuple(W[:,k]) for k in range(D)]
    if any(c[0]==0 for c in cols): continue
    if len({c for c in cols}|{tuple(-x for x in c) for c in cols})<2*D: continue
    # roots: y=w1, x=w2/(2y)
    roots=[sp.Rational(c[1])/(2*c[0])+sp.I*c[0] for c in cols]
    Q=sp.expand(sp.prod([t-r for r in roots]))
    cnt=0
    for F in itertools.product([0,1],repeat=D):
        z=sp.expand(sp.prod([(t-sp.conjugate(roots[k])) if F[k] else (t-roots[k]) for k in range(D)]))
        dz=sp.expand(z-Q)
        if dz==0 or sp.degree(dz,t)<=K: cnt+=1
    print("exact cluster size through Q with deg<=K=%d:"%K,cnt,"roots",roots); break
