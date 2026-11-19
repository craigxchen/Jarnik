# h^0(B, jM), M=-K_B = 5h - E3 - 2E4 - 3E5 on the wonderful model of the braid arrangement
# (6 labels 0..5 mod translation, P^4 with coords y_1..y_5, y_0=0).
# Section = form of degree 5j with mult >= j (3-clusters), 2j (4-clusters), 3j (5-clusters).
# Conditions are imposed EXACTLY (coordinate change to cluster-adapted coords), rank mod p.
# Also variant "V_j": only 3- and 4-cluster conditions (Astra's induction space).
import sys, itertools
from math import comb
import flint
p = int(sys.argv[2]) if len(sys.argv)>2 else 1000003
j = int(sys.argv[1])
d = 5*j
def monos(d,n=5):
    if n==1: yield (d,); return
    for a in range(d,-1,-1):
        for r in monos(d-a,n-1): yield (a,)+r
M = list(monos(d)); idx = {m:i for i,m in enumerate(M)}; N=len(M)
labels = range(6)
# transform for base a: new coords z_c = y_c - y_a, c != a (z_0 = -y_a).
# old y_c (=y_c - y_0) = z_c - z_0 for c not in {0,a}; y_a = -z_0.
# new variable order: list of c != a in increasing order (includes 0 when a!=0).
def transform_rows(a, needed):
    """return dict newmono -> {oldindex: coeff} for new monomials in 'needed'"""
    newvars = [c for c in labels if c != a]  # length 5
    pos = {c:i for i,c in enumerate(newvars)}
    out = {}
    for oi,m in enumerate(M):
        # m exponents for y_1..y_5
        if a == 0:
            nm = tuple(m[pos[c]-0] if False else m[c-1] for c in newvars)
            if nm in needed: out.setdefault(nm,{})[oi]=1
            continue
        # factors: for c in 1..5, c!=a: (z_c - z_0)^{m[c-1]}; for c=a: (-z_0)^{m[a-1]}
        others = [c for c in range(1,6) if c != a]
        ranges = [range(m[c-1]+1) for c in others]
        for ks in itertools.product(*ranges):
            coef = (-1)**m[a-1]
            e0 = m[a-1]
            nm = [0]*5
            for c,k in zip(others,ks):
                coef *= comb(m[c-1],k)*(-1)**(m[c-1]-k)
                nm[pos[c]] = k; e0 += m[c-1]-k
            nm[pos[0]] = e0
            nm = tuple(nm)
            if nm in needed:
                dd = out.setdefault(nm,{}); dd[oi] = dd.get(oi,0)+coef
    return out, newvars
def build(with_points=True):
    rows = []
    sizes = [(3,j),(4,2*j)] + ([(5,3*j)] if with_points else [])
    for a in labels:
        newvars = [c for c in labels if c != a]
        pos = {c:i for i,c in enumerate(newvars)}
        needed = set(); clusters=[]
        for s,mult in sizes:
            for S in itertools.combinations(labels,s):
                if min(S)!=a: continue
                Sp = [pos[c] for c in S if c!=a]
                clusters.append(S)
                for nm in monos(d):
                    if sum(nm[i] for i in Sp) < mult: needed.add(nm)
        if not needed: continue
        out,_ = transform_rows(a, needed)
        for nm in needed:
            r = out.get(nm,{})
            if r: rows.append(r)
    return rows
def rank_kernel(rows):
    A = flint.nmod_mat(len(rows), N, p)
    for i,r in enumerate(rows):
        for k,v in r.items(): A[i,k] = v % p
    rk = A.rank()
    return len(rows), rk, N-rk
for wp in (True, False):  # use for j<=3; for j>=4 use h0M3.py (column elimination, memory)
    rows = build(wp)
    nr, rk, ker = rank_kernel(rows)
    print(f"j={j} deg={d} monomials={N} pointconds={wp} rows={nr} rank={rk} kernel={ker}", flush=True)
