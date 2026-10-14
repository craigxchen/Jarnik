# Syl(16,2) plane image P (7 columns of a 3-dim subspace V), in its own ring Q[y_1..y_7]:
# J(P) (mod p, upper bound) vs J_w(P) + line pullbacks + families of <=b directions from natural candidates.
import numpy as np, itertools, sys, time, math
from flint import nmod_mat
from rengc import REngC, PR
from rlib import J_exact, monos
from syl16_base import COLS, pts
DT=int(sys.argv[1]); bmax=int(sys.argv[2])
A=np.load("syl16_A.npy").astype(np.int64)
Vi=[pts.index(p) for p in pts if p[0]==0]; Va=[pts[i] for i in Vi]
P=np.unique(A@COLS[Vi].T.astype(np.int64),axis=0); print("|P|",len(P))
E=REngC(7,DT,cap=None)
Zm=[]  # evaluation for J(P) mod p
blocks=[]
for D in range(DT+1):
    M=np.array(E.ms[D],dtype=np.int64)
    val=np.ones((len(P),len(M)),dtype=object)
    for j in range(7):
        val=val*np.power(P[:,j:j+1].astype(object),M[None,:,j])
    blocks.append(np.array(val%PR,dtype=np.int64))
JP=[];prev=0
for D in range(DT+1):
    S=np.concatenate(blocks[:D+1],axis=1)
    r=nmod_mat(S.shape[0],S.shape[1],[int(x) for x in S.flatten()],PR).rank()
    JP.append(len(E.ms[D])-(r-prev)); prev=r
print("dim J(P)_D mod p (upper bounds)",JP,flush=True)
class Acc:
    def __init__(s): s.B=None; s.r=0
    def add(s,rows):
        if not rows: return 0
        M=np.vstack(([s.B] if s.B is not None else [])+[np.array(rows)])
        nm=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR); R,rk=nm.rref()
        k=rk-s.r
        if k>0: s.B=np.array([[int(R[i,j]) for j in range(M.shape[1])] for i in range(rk)],dtype=np.int64); s.r=rk
        return k
def fam_rows(U,Dt,skip_grid=True):
    U=[np.array(u,dtype=np.int64) for u in U]
    img=np.unique(np.stack([P@u for u in U],axis=1),axis=0)
    sizes=[len(np.unique(img[:,i])) for i in range(len(U))]
    if skip_grid and np.prod(sizes)==len(img): return None,len(img)
    J=J_exact([tuple(int(x) for x in p) for p in img],Dt)
    rows=[]
    for D in range(Dt+1):
        tm=monos(len(U),D)
        for vec in J[D]:
            g=E.pullback(vec,tm,U)
            if g.any(): rows+=E.ideal_rows(g,D,Dt)
    return rows,len(img)
# candidate directions in y-coordinates (integer vectors, scaled)
cand={}
for i in range(7):
    e=[0]*7; e[i]=1; cand['y%d'%i]=e
lines=[]
for i,j in itertools.combinations(range(7),2):
    c=tuple((a+b)%2 for a,b in zip(Va[i],Va[j])); k=Va.index(c); l=tuple(sorted((i,j,k)))
    if l not in lines: lines.append(l)
for l in lines:
    e=[0]*7
    for i in l: e[i]=1
    cand['z%s'%(l,)]=e      # = 2 z_l
cand['zV']=[1]*7           # = 4 z_V
# coset sums kappa_c: c ranges over nonzero classes ~ nonzero elements of dual; kappa_c*4 = sum_{a.c=1} y_a - sum_{a.c=0} y_a
for c in Va:
    cand['k%s'%(c,)]=[1 if sum(x*y for x,y in zip(a,c))%2 else -1 for a in Va]
names=list(cand); print("candidates",len(names))
# width part: all primitive integer dirs in [-1,1]^7 of the y-ring, plus candidates
acc=Acc()
Wd=[v for v in itertools.product([-1,0,1],repeat=7) if any(v)]+[tuple(v) for v in cand.values()]
for v in Wd:
    s=len(np.unique(P@np.array(v)))
    if s<=DT: acc.add(E.ideal_rows(E.pullback([1],[(s,)],[list(v)]),s,DT))
print("D",DT,"J_w(P)",acc.r,"J(P)",JP[DT],flush=True)
for l in lines:
    U=[cand['y%d'%i] for i in l]; rows,sz=fam_rows(U,DT); acc.add(rows)
print("+lines",acc.r,flush=True)
for b in range(2,bmax+1):
    t0=time.time(); n=0
    for F in itertools.combinations(names,b):
        U=[cand[f] for f in F]
        if np.linalg.matrix_rank(np.array(U))<b: continue
        rows,sz=fam_rows(U,DT)
        if rows is None: continue
        n+=1
        k=acc.add(rows)
        if k: print(" b",b,F,"|img|",sz,"new",k,"total",acc.r,flush=True)
    print("after all b=%d families (%d nongrid):"%(b,n),acc.r,"time",time.time()-t0,flush=True)
