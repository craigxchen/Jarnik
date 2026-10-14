# Syl(16,2), Singer-invariant part of degree-7 ideals (mod p, evidence): author's width set + lines, then + affine 4-sums.
import numpy as np, itertools, time, sys
from flint import nmod_mat
from syl16_base import COLS, pts
from reng import REng, PR
from rlib import J_exact, monos
A=np.load("syl16_A.npy").astype(np.int64)
d=15; DT=7
def mult_t(p):  # p=(c0,c1,c2,c3) coefficients of 1,t,t^2,t^3; multiply by t mod t^4+t+1
    c0,c1,c2,c3=p
    return (c3, (c0+c3)%2, c1, c2)
shift=[pts.index(mult_t(p)) for p in pts]   # coordinate j -> shift[j]
assert sorted(shift)==list(range(15))
t0=time.time(); E=REng(d,DT,shift); print("engine",time.time()-t0, "orbits deg7",E.orb[DT][1],flush=True)
class Acc:
    def __init__(s): s.B=None; s.r=0
    def add(s,rows):
        if not rows: return 0
        M=np.vstack(([s.B] if s.B is not None else [])+[np.array(rows)])
        nm=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR); R,rk=nm.rref()
        k=rk-s.r
        if k>0: s.B=np.array([[int(R[i,j]) for j in range(M.shape[1])] for i in range(rk)],dtype=np.int64); s.r=rk
        return k
def P0rows(g,D): return [E.P0(v,DT) for v in E.ideal_rows(g,D,DT)]
def orbit_reps(vecs):
    reps={}
    for v in vecs:
        key=tuple(v); best=key; cur=key
        for t in range(15):
            cur=tuple(cur[shift.index(j)] for j in range(15)); best=min(best,cur)
        reps.setdefault(best,np.array(v))
    return list(reps.values())
def width_gen(v):
    s=len(np.unique(A@v)); return E.pullback([1],[(s,)],[v]),s
acc=Acc()
one=np.ones(15,dtype=np.int64)
# columns
for v in orbit_reps([COLS[i].astype(np.int64) for i in range(15)]):
    g,s=width_gen(v); acc.add(P0rows(g,s))
print("cols",acc.r,flush=True)
g,s=width_gen(one); acc.add(P0rows(g,s)); print("+sigma^%d"%s,acc.r,flush=True)
lines=[]
for i,j in itertools.combinations(range(15),2):
    c=tuple((a+b)%2 for a,b in zip(pts[i],pts[j])); k=pts.index(c); l=tuple(sorted((i,j,k)))
    if l not in lines: lines.append(l)
half=[(COLS[a]+COLS[b]+COLS[c]).astype(np.int64)//2 for a,b,c in lines]
for v in orbit_reps(half):
    g,s=width_gen(v); acc.add(P0rows(g,s))
print("+half-sums",acc.r,flush=True)
planes=[]
for a in pts:
    V=[i for i,p in enumerate(pts) if sum(x*y for x,y in zip(a,p))%2==0]
    planes.append(V)
quarter=[COLS[V].sum(0).astype(np.int64)//4 for V in planes]
for v in orbit_reps(quarter):
    g,s=width_gen(v); acc.add(P0rows(g,s))
print("+quarter-sums (author's J_w)",acc.r,flush=True)
# line families: orbit reps of lines
lreps={}
for l in lines:
    key=tuple(sorted(l)); best=key; cur=key
    for t in range(15):
        cur=tuple(sorted(shift[i] for i in cur)); best=min(best,cur)
    lreps.setdefault(best,l)
for l in lreps.values():
    U=[COLS[i].astype(np.int64) for i in l]
    img=np.unique(np.stack([A@u for u in U],axis=1),axis=0)
    J=J_exact([tuple(int(x) for x in p) for p in img],DT)
    for D in range(DT+1):
        for vec in J[D]:
            g=E.pullback(vec,monos(3,D),U)
            if g.any(): acc.add(P0rows(g,D))
print("+lines",acc.r,flush=True)
base=acc.r
Bbase=acc.B.copy()
W7=np.load("syl16_w7dirs.npy")
aff=[v for v in W7 if np.abs(v).sum()==4]
print("support-4 width-7 dirs",len(aff), "example", aff[0].tolist(), [pts[i] for i in np.nonzero(aff[0])[0]])
for v in orbit_reps(aff):
    g,s=width_gen(v.astype(np.int64)); k=acc.add(P0rows(g,s)); print("  affine 4-sum",v.tolist(),"width",s,"new",k,"total",acc.r,flush=True)
print("author J_w+lines:",base," + affine 4-sums:",acc.r,flush=True)
import pickle
Baff=acc.B.copy(); raff=acc.r
# plane generators mod p (evidence): J(P) via mod-p nullspace of evaluation
V0=planes[0]; U=[COLS[i].astype(np.int64) for i in V0]
img=np.unique(np.stack([A@u for u in U],axis=1),axis=0); print("plane image",len(img),flush=True)
allm=[]; 
for D in range(DT+1): allm+=[(D,m) for m in monos(7,D)]
def ev(p,m):
    r=1
    for a,b in zip(p,m): r*=int(a)**b
    return r%PR
gens=[]
for D in range(5,DT+1):
    top=[m for DD,m in allm if DD==D]; low=[m for DD,m in allm if DD<D]
    M=nmod_mat([[ev(p,m) for m in top+low] for p in img],PR)
    N=M.nullspace()[0]; nul=M.nullspace()[1]
    vecs=[[int(N[i,j]) for i in range(len(top))] for j in range(nul)]
    T=nmod_mat(vecs,PR) if vecs else None
    if T is None: continue
    R,rk=T.rref()
    for i in range(rk):
        gens.append((D,[int(R[i,j]) for j in range(len(top))],top))
    print(" plane J_%d mod p dim"%D,rk,flush=True)
PR_rows=[]
for D,vec,top in gens:
    g=E.pullback(vec,top,U)
    if g.any(): PR_rows+=P0rows(g,D)
print("plane P0 rows",len(PR_rows),flush=True)
def addplane(acc):
    acc.add(PR_rows)
acc2=Acc(); acc2.B=Bbase.copy(); acc2.r=base
addplane(acc2); print("J_w+lines+planes:",acc2.r,flush=True)
acc.B=Baff; acc.r=raff; addplane(acc); print("J_w+lines+affine+planes:",acc.r,flush=True)
