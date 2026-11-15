# Exhaustive: base-point codes C (0 in C) of real rank r<=4 in {0,1}^D, distances >= L+1,
# with balance: every difference has >= ceil(L/2) entries of each sign after a sign pattern sigma.
# Info-set reduction: C ~ S subset {0,1}^r spanning, columns = linear functionals l with l(S) in {0,1}.
import itertools, numpy as np, sys
from scipy.optimize import milp, LinearConstraint, Bounds
def functionals(S, r):
    S=np.array(S); # find basis inside S
    import itertools
    for B in itertools.combinations(range(len(S)), r):
        Mb=S[list(B)]
        if abs(np.linalg.det(Mb))>0.5: break
    inv=np.linalg.inv(Mb)
    out=[]
    for vals in itertools.product([0,1],repeat=r):
        l=inv@np.array(vals,float)
        v=S@l
        if np.all(np.abs(v-np.round(v))<1e-9) and np.all((np.round(v)==0)|(np.round(v)==1)):
            out.append(np.round(v).astype(int))
    return out  # columns as value vectors on S
def feasible(S, r, D, L, balance=True, want=None):
    cols=functionals(S,r)
    n=len(S); pairs=list(itertools.combinations(range(n),2))
    # variables: for each col, two signs
    nv=2*len(cols)
    A=[];lb=[];ub=[]
    A.append([1]*nv); lb.append(D); ub.append(D)
    hb=(L+1)//2  # ceil(L/2)
    for (x,y) in pairs:
        dist=[];pos=[];neg=[]
        for c in cols:
            dd=c[x]-c[y]
            for s in (1,-1):
                dist.append(abs(dd)); pos.append(1 if s*dd==1 else 0); neg.append(1 if s*dd==-1 else 0)
        A.append(dist); lb.append(L+1); ub.append(np.inf)
        if balance:
            A.append(pos); lb.append(hb); ub.append(np.inf)
            A.append(neg); lb.append(hb); ub.append(np.inf)
    # j=1 Gordan-type necessary condition is checked afterwards
    res=milp(c=np.zeros(nv),constraints=LinearConstraint(np.array(A,float),lb,ub),
             integrality=np.ones(nv),bounds=Bounds(0,D))
    if res.status==0:
        m=np.round(res.x).astype(int)
        return [(cols[i//2],1 if i%2==0 else -1,m[i]) for i in range(nv) if m[i]>0]
    return None
def run(D,L,r=4,balance=True,smin=None):
    vecs=[v for v in itertools.product([0,1],repeat=r) if any(v)]
    best=0;ex=None;count={}
    for s in range(2**r, r, -1):
        if smin and s<smin: break
        found=0
        for sub in itertools.combinations(vecs,s-1):
            S=[(0,)*r]+list(sub)
            if np.linalg.matrix_rank(np.array(S))<r: continue
            f=feasible(S,r,D,L,balance)
            if f is not None:
                found+=1
                if ex is None: ex=(S,f)
        count[s]=found
        if found: return s,count,ex
    return None
if __name__=="__main__":
    D=int(sys.argv[1]); L=int(sys.argv[2]); bal=sys.argv[3]=='1'; r=int(sys.argv[4]) if len(sys.argv)>4 else 4
    out=run(D,L,r,bal)
    print("D",D,"L",L,"balance",bal,"rank",r,"max |C|",out[0],"num S (labelled) at max",out[1][out[0]])
    S,f=out[2]; print("example S",S); 
    for c,s,m in f: print("  col",list(c),"sigma",s,"mult",m)
