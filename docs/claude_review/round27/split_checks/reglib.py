# Library: regularity index reg(S) of a finite point set S in M = {m in Z^5 : sum m = 0} (coords m1..m4),
# = max order at e of a Laurent polynomial supported on S; toric pole orders pi_B(S) = max_{m in S} (-<m,1_B>).
import itertools, numpy as np
from flint import nmod_mat
n=5
Bs=[frozenset(c) for r in range(1,n) for c in itertools.combinations(range(n),r)]
idx={B:i for i,B in enumerate(Bs)}
comp=lambda B: frozenset(range(n))-B
Bmat=np.array([[1 if a in B else 0 for a in range(5)] for B in Bs])  # 30 x 5
def full(m4): return tuple(list(m4)+[-sum(m4)])
def pis(S):
    S5=np.array([full(m) for m in S]); return (-(S5@Bmat.T)).max(axis=0)   # length 30
_monos={}
def monos(d):
    if d not in _monos:
        _monos[d]=[a for t in range(d+1) for a in itertools.product(range(t+1),repeat=4) if sum(a)==t]
        _monos[d]=sorted(set(_monos[d]),key=lambda a:(sum(a),a))
    return _monos[d]
def reg(S,p=2**31-1,dmax=60):
    """min d with rank of [m^alpha]_{m in S, |alpha|<=d} = |S| (mod p). Upper bound for true reg when mod-p rank full."""
    S=list(S); N=len(S)
    lo=0
    for d in range(0,dmax+1):
        ms=monos(d)
        if len(ms)<N: continue
        rows=[[pow(int(m[0])%p,a[0],p)*pow(int(m[1])%p,a[1],p)%p*pow(int(m[2])%p,a[2],p)%p*pow(int(m[3])%p,a[3],p)%p for a in ms] for m in S]
        M=nmod_mat(rows,p)
        if M.rank()==N: return d
    return None
def lattice_points(p_vec, box=None):
    """points m in Z^4 with <m,1_B> >= -p_B for all B with p_B not None"""
    from scipy.optimize import linprog
    A=[];b=[]
    for i,B in enumerate(Bs):
        if p_vec[i] is None: continue
        coef=np.zeros(4)
        for a in B:
            if a<4: coef[a]+=1
            else: coef-=1
        A.append(-coef); b.append(p_vec[i])
    A=np.array(A,float); b=np.array(b,float)
    lo=[];hi=[]
    for i in range(4):
        c=np.zeros(4); c[i]=1
        r1=linprog(c,A_ub=A,b_ub=b,bounds=[(None,None)]*4,method='highs')
        r2=linprog(-c,A_ub=A,b_ub=b,bounds=[(None,None)]*4,method='highs')
        if r1.status!=0 or r2.status!=0: return None
        lo.append(int(np.floor(r1.fun+1e-9))); hi.append(int(np.ceil(-r2.fun-1e-9)))
    pts=[]
    for m in itertools.product(*[range(lo[i],hi[i]+1) for i in range(4)]):
        if np.all(A@np.array(m)<=b+1e-9): pts.append(m)
    return pts
