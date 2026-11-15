# independent check of Prop 3.3: codes C = V cap {0,1}^D, dim V <= r, 0 in C, min distance >= theta*D
import itertools, numpy as np, sys
from fractions import Fraction
from scipy.optimize import linprog
r=int(sys.argv[1])
vecs=[v for v in itertools.product([0,1],repeat=r) if any(v)]
def rank(M): return np.linalg.matrix_rank(np.array(M,dtype=float)) if M else 0
best={}  # size -> max t*
maxsize_at={}
def functionals(S):
    # find basis in S
    B=[];
    for v in S:
        if rank(B+[v])>len(B): B.append(v)
        if len(B)==r: break
    if len(B)<r: return None
    Binv=np.linalg.inv(np.array(B,dtype=float))
    out=set()
    for vals in itertools.product([0,1],repeat=r):
        if not any(vals): continue
        l=Binv@np.array(vals,dtype=float)  # l such that B l = vals
        ok=True; img=[]
        for s in S:
            x=float(np.dot(s,l))
            if abs(x)<1e-9: img.append(0)
            elif abs(x-1)<1e-9: img.append(1)
            else: ok=False;break
        if ok: out.add(tuple(img))
    return list(out)
res=[]
for mask in range(1,2**len(vecs)):
    S=[(0,)*r]+[vecs[i] for i in range(len(vecs)) if mask>>i&1]
    if len(S)<r+1: continue
    F=functionals(S)
    if F is None: continue
    n=len(S); pairs=list(itertools.combinations(range(n),2))
    k=len(F)
    # variables w_1..w_k, t ; maximize t
    A=[];b=[]
    for (i,j) in pairs:
        row=[-abs(f[i]-f[j]) for f in F]+[1.0]
        A.append(row); b.append(0.0)
    Aeq=[[1.0]*k+[0.0]]; beq=[1.0]
    c=[0.0]*k+[-1.0]
    sol=linprog(c,A_ub=A,b_ub=b,A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*k+[(None,None)],method='highs')
    t=-sol.fun
    best[n]=max(best.get(n,0),t)
for n in sorted(best): print(r,n,round(best[n],6))
