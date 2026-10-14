# referee's own library for M8 (independent of refJ_checks/m8lib.py)
import itertools, math, random
from flint import fmpz_mat, nmod_mat
P = 2**61-1  # large prime for mod-p lower bounds
NV = 7
def dot(a,b): return bin(a&b).count('1')&1
# points x=1..7 (nonzero of F_2^3); T_a = {x: a.x=1}
T = {a:[x for x in range(1,8) if dot(a,x)] for a in range(1,8)}
def col(a): return tuple(1 if x in T[a] else 0 for x in range(1,8))   # coordinate i <-> x=i+1
def build_A(L=2):
    A=[]
    for k in itertools.product((-1,0,1),repeat=7):
        if all(abs(sum(k[x-1] for x in T[a]))<=L for a in range(1,8)): A.append(k)
    return A
def monos(nv,e):
    out=[]
    for c in itertools.combinations_with_replacement(range(nv),e):
        ex=[0]*nv
        for i in c: ex[i]+=1
        out.append(tuple(ex))
    return out
def monos_le(nv,e): return [m for t in range(e+1) for m in monos(nv,t)]
def evmono(m,p): 
    r=1
    for e_,y_ in zip(m,p):
        if e_: r*=y_**e_
    return r
def rank_exact(rows):
    if not rows: return 0
    return fmpz_mat(rows).rank()
def rank_modp(rows,p=P):
    if not rows: return 0
    return nmod_mat([[x%p for x in r] for r in rows],p).rank()
# polynomials as dict exponent-tuple -> int
def pmul(f,g):
    h={}
    for a,x in f.items():
        for b,y in g.items():
            c=tuple(i+j for i,j in zip(a,b)); h[c]=h.get(c,0)+x*y
    return {c:v for c,v in h.items() if v}
def ppow(f,t):
    r={tuple([0]*len(next(iter(f)))):1}
    for _ in range(t): r=pmul(r,f)
    return r
def linform(v):
    n=len(v); return {tuple(1 if i==j else 0 for i in range(n)):int(v[j]) for j in range(n) if v[j]}
def compose(f, V):
    """f: poly in b vars (dict); V: list of b integer vectors in Z^7; returns f(<x,V_1>,...,<x,V_b>) homogeneous parts as dict"""
    lfs=[linform(v) for v in V]
    out={}
    cache={}
    for ex,cf in f.items():
        term={tuple([0]*NV):cf}
        for i,e in enumerate(ex):
            if e:
                key=(i,e)
                if key not in cache: cache[key]=ppow(lfs[i],e)
                term=pmul(term,cache[key])
        for c,v in term.items(): out[c]=out.get(c,0)+v
    return {c:v for c,v in out.items() if v}
def vec(f, basis_index):
    v=[0]*len(basis_index)
    for c,x in f.items(): v[basis_index[c]]+=x
    return v
def kernel_q(M, ncols):
    """exact kernel basis (list of lists of fmpq) of rational matrix M (list of rows)"""
    from flint import fmpq_mat, fmpq
    if not M: return [[fmpq(1) if i==j else fmpq(0) for i in range(ncols)] for j in range(ncols)]
    R,rk=fmpq_mat(M).rref()
    piv=[]; 
    for i in range(rk):
        for t in range(ncols):
            if R[i,t]!=0: piv.append(t); break
    free=[t for t in range(ncols) if t not in piv]
    K=[]
    for f in free:
        v=[fmpq(0)]*ncols; v[f]=fmpq(1)
        for i,pc in enumerate(piv): v[pc]=-R[i,f]
        K.append(v)
    return K
def leading_ideal_gens(Y, nv, maxdeg):
    """exact basis of J(Y)_e, e<=maxdeg, as integer dict polys. J(Y)_e={f in R_e: ev f in ev R_<e}"""
    from flint import fmpq_mat
    res={}
    for e in range(maxdeg+1):
        Me=monos(nv,e); Ml=monos_le(nv,e-1) if e>0 else []
        cols=Me+Ml
        M=[[evmono(m,y) for m in cols] for y in Y]
        K=kernel_q(M,len(cols))
        vecs=[v[:len(Me)] for v in K]
        gens=[]
        if vecs:
            Qr,rk=fmpq_mat(vecs).rref()
            for i in range(rk):
                row=[Qr[i,t] for t in range(len(Me))]
                den=1
                for x in row: den=den*int(x.q)//math.gcd(den,int(x.q))
                gens.append({Me[t]:int(row[t]*den) for t in range(len(Me)) if row[t]!=0})
        res[e]=gens
    return res
def rho(Y,nv):
    """interpolation degree: smallest e with rank ev(R_<=e) = |Y| (exact)"""
    e=0
    while True:
        r=rank_exact([[evmono(m,y) for m in monos_le(nv,e)] for y in Y])
        if r==len(Y): return e
        e+=1
