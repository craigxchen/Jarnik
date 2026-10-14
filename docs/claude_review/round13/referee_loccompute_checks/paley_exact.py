# Exact (fmpz) rank of P_0(L_7) for Paley-12 L=2, for
#   L1 = J_w(R=1) + all 5-column families           (author: 82)
#   L2 = L1 + the sigma+3-coordinate families         (referee: >= 87 mod p)
# Completeness of the 5-column generators: congruence families -> exact J(P) (fmpz) in degrees <= 7;
# near-full families -> mod-p check that J(P)_{<=7} = (y_q^5)_{<=7} (mod-p rank of evaluation gives an upper
# bound on dim_Q J(P), the ideal (y_q^5) is inside J(P) over Q, equal dims => equal over Q).
import numpy as np, itertools, time, sys
from math import comb
from flint import fmpz_mat, nmod_mat
from paley_base import build, COLS
from rlib import J_exact, monos
P=2147483629
A=build(2); d=11; DT=7
# ---- exact polynomial arithmetic in R (dict exponent->int), exponent cap 2
def mul(f,g,Dmax=DT):
    out={}
    for e1,c1 in f.items():
        for e2,c2 in g.items():
            e=tuple(a+b for a,b in zip(e1,e2))
            if max(e)>2: continue
            out[e]=out.get(e,0)+c1*c2
    return {e:c for e,c in out.items() if c}
def lin(v): return {tuple(1 if j==i else 0 for j in range(d)):int(v[i]) for i in range(d) if v[i]}
def pull(vec,tm,U):
    Ls=[lin(u) for u in U]; cache={}
    def pw(i,k):
        if (i,k) in cache: return cache[(i,k)]
        r={tuple([0]*d):1} if k==0 else mul(pw(i,k-1),Ls[i]); cache[(i,k)]=r; return r
    out={}
    for c,e in zip(vec,tm):
        if not c: continue
        t={tuple([0]*d):int(c)}
        for i,k in enumerate(e):
            if k: t=mul(t,pw(i,k))
        for m,cc in t.items(): out[m]=out.get(m,0)+cc
    return {m:c for m,c in out.items() if c}
ms7=monos(d,DT,2)
def canon(m): return min(tuple(m[(j+t)%d] for j in range(d)) for t in range(d))
orbits=sorted({canon(m) for m in ms7}); oidx={o:i for i,o in enumerate(orbits)}
print("degree-7 orbits",len(orbits))
def P0row(f):
    r=[0]*len(orbits)
    for m,c in f.items(): r[oidx[canon(m)]]+=c
    return r
def rows_of(g,D):
    out=[]
    for mu in monos(d,DT-D,2):
        f=mul(g,{mu:1})
        if f: out.append(P0row(f))
    return out
def rank_exact(rows):
    return fmpz_mat(rows).rank() if rows else 0
rows=[]
# J_w(R=1): width<=7 directions, orbit reps (from paley_base scan)
V=np.load("paley_V.npy"); w=np.load("paley_w.npy")
def canv(v): return max(tuple(int(x) for x in s*np.roll(v,t)) for t in range(11) for s in (1,-1))
reps={}
for v in V[w<=7]: reps.setdefault(canv(v),v)
for v in reps.values():
    s=int(len(np.unique(A@v))); g=pull([1],[(s,)],[v])
    if g: rows+=rows_of(g,s)
print("J_w rows",len(rows),"exact rank",rank_exact(rows),flush=True)
# 5-column families
seen=set(); t0=time.time()
for F in itertools.combinations(range(11),5):
    c=min(tuple(sorted((i+t)%11 for i in F)) for t in range(11))
    if c in seen: continue
    seen.add(c)
    U=[COLS[i] for i in c]
    img=np.unique(np.stack([A@u for u in U],axis=1),axis=0)
    if len(img)>=3000:
        # near-full: mod-p check J(P)_{<=7} == (y^5)_{<=7}
        prev=0; ok=True
        for D in range(DT+1):
            M=[]
            ex=[m for DD in range(D+1) for m in monos(5,DD)]
            Ev=nmod_mat([[int(np.prod([int(p[j])**m[j] for j in range(5)]))%P for m in ex] for p in img],P)
            r=Ev.rank(); H=r-prev; prev=r
            jd=comb(D+4,4)-H
            want=0 if D<5 else 5*comb(D-5+4,4)
            if jd!=want: ok=False
        print("near-full",c,len(img),"J(P)_{<=7}=(y^5):",ok,flush=True)
        if not ok: sys.exit("completeness failed")
        continue
    J=J_exact([tuple(int(x) for x in p) for p in img],DT)
    for D in range(DT+1):
        for vec in J[D]:
            g=pull(vec,monos(5,D),U)
            if g: rows+=rows_of(g,D)
    print("congruence",c,len(img),"J dims",[len(J[D]) for D in range(DT+1)],"rows",len(rows),"exact rank",rank_exact(rows),time.time()-t0,flush=True)
r1=rank_exact(rows); print("EXACT dim_Q P0(L_col(5))_7 =",r1,flush=True)
one=np.ones(11,dtype=np.int64)
for S in [(0,1,5),(0,2,6),(0,1,2),(0,1,3)]:
    U=[one]+[np.array([1 if j==i else 0 for j in range(11)]) for i in S]
    img=np.unique(np.stack([A@u for u in U],axis=1),axis=0)
    J=J_exact([tuple(int(x) for x in p) for p in img],DT)
    for D in range(DT+1):
        for vec in J[D]:
            g=pull(vec,monos(4,D),U)
            if g: rows+=rows_of(g,D)
    print("after sigma+coords",S,"exact rank",rank_exact(rows),flush=True)
