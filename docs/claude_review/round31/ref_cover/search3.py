# Evidence search: windows of consecutive lattice points on x^2+y^2=N, all 5-subsets, gcd-reduce,
# compute index I of valuation-column lattice and g_eff = log H / (-log Delta). Small g_eff = good approximation.
import sys, math, cmath
from itertools import combinations, product
from math import gcd, log
sys.path.insert(0,'.')
from gauss import *
from lattice import lattice_index

PR=[5,13,17,29,37,41,53,61,73]
def reps(fac):
    # all z with |z|^2 = prod p^e  (all split), primitive or not
    pts=[(1,0)]
    for p,e in fac:
        pi=PI[p]; new=[]
        for a in range(e+1):
            f=gmul(gpow(pi,a),gpow(gconj(pi),e-a))
            for w in pts: new.append(gmul(w,f))
        pts=new
    out=[]
    for w in pts:
        for u in UNITS: out.append(gmul(w,u))
    return out

def index_of(cols):
    # cols: list of 4-vectors; return (rank, index)
    import sympy as sp
    M=sp.Matrix(cols).T
    r=M.rank()
    if r<4: return r,None
    I=0
    for c in combinations(range(len(cols)),4):
        I=gcd(I,int(abs(M[:,list(c)].det())))
    return 4,I

def gen_N(limit):
    res=[]
    def rec(i,N,fac):
        if i==len(PR):
            if len(fac)>=1: res.append((N,list(fac)))
            return
        p=PR[i]; e=0; M=N
        while M<=limit:
            rec(i+1,M,fac+([(p,e)] if e else []))
            M*=p; e+=1
    rec(0,1,[])
    return res
limit=float(sys.argv[1]); GMAX=float(sys.argv[2]); HMIN=int(sys.argv[3])
Ns=gen_N(limit)
Ns=[(N,f) for N,f in Ns if math.prod(e+1 for _,e in f)>=96]
print("N count",len(Ns),file=sys.stderr)
best={}  # key: index class -> list of (geff, N, H')
WIN=8
for N,fac in Ns:
    pts=reps(fac)
    pts=sorted(set(pts), key=lambda w: math.atan2(w[1],w[0]))
    n=len(pts); R=math.sqrt(N)
    seen=set()
    for s in range(n):
        win=[pts[(s+t)%n] for t in range(min(WIN,n))]
        for S in combinations(win,5):
            key=frozenset(S)
            if key in seen: continue
            seen.add(key)
            # diameter
            D=max(abs(complex(*a)-complex(*b)) for a,b in combinations(S,2))
            g=S[0]
            for w in S[1:]: g=ggcd(g,w)
            ng=gnorm(g)
            Hp=N//ng
            if Hp<HMIN: continue
            Delta=D/R
            geff=log(Hp)/(-log(Delta))
            if geff>GMAX: continue
            red=[gdivmod_exact(w,g) for w in S]
            cols=[]
            for p,e in fac:
                if Hp%p: continue
                pi=PI[p]
                a=[val(w,pi) for w in red]
                cols.append([a[j]-a[4] for j in range(4)])
            if not cols: continue
            r,I=lattice_index(cols)
            cls = ('rank%d'%r) if I is None else ('I=%d'%I)
            best.setdefault(cls,[]).append((round(geff,4),N,Hp))
for cls in sorted(best):
    L=sorted(set(best[cls]))[:4]
    print(cls, len(best[cls]), L)
