import numpy as np, itertools, sys, time
from flint import nmod_mat
from reng import REng, PR
from rlib import J_exact, monos
L=int(sys.argv[1]); DT=int(sys.argv[2])
pts=[p for p in itertools.product([0,1],repeat=3) if any(p)]
d=7
COLS=np.array([[1 if sum(a*x for a,x in zip(av,p))%2 else 0 for p in pts] for av in pts],dtype=np.int64)
K=np.array(list(itertools.product([-1,0,1],repeat=d)),dtype=np.int64)
A=K[np.all(np.abs(K@COLS.T)<=L,axis=1)]
print("|A|",len(A))
E=REng(d,DT)
# J(A) dims mod p (upper bounds)
Zr=(A==0).astype(np.int64); Nr=(A==-1).astype(np.int64)
blocks=[]
for D in range(DT+1):
    M=np.array(E.ms[D],dtype=np.int64)
    zero=(Zr@(M>0).T.astype(np.int64))>0; par=(Nr@(M%2==1).T.astype(np.int64))%2
    blocks.append(np.where(zero,0,1-2*par))
JA=[];prev=0
for D in range(DT+1):
    S=np.concatenate(blocks[:D+1],axis=1)
    r=nmod_mat(S.shape[0],S.shape[1],[int(x)%PR for x in S.flatten()],PR).rank()
    JA.append(len(E.ms[D])-(r-prev)); prev=r
print("dim J(A)_D (mod p upper bounds)",JA)
class Acc:
    def __init__(s,D): s.D=D; s.B=None; s.r=0
    def add(s,rows):
        if not rows: return 0
        M=np.vstack(([s.B] if s.B is not None else [])+[np.array(rows)])
        nm=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR); R,rk=nm.rref()
        k=rk-s.r
        if k>0: s.B=np.array([[int(R[i,j]) for j in range(M.shape[1])] for i in range(rk)],dtype=np.int64); s.r=rk
        return k
def fam_gens(U):
    U=[np.array(u,dtype=np.int64) for u in U]
    img=np.unique(np.stack([A@u for u in U],axis=1),axis=0)
    sizes=[len(np.unique(img[:,i])) for i in range(len(U))]
    grid=np.prod(sizes)==len(img)
    J=J_exact([tuple(int(x) for x in p) for p in img],DT)
    gens=[]
    for D in range(DT+1):
        tm=monos(len(U),D)
        for vec in J[D]:
            g=E.pullback(vec,tm,U)
            if g.any(): gens.append((g,D))
    return gens,len(img),grid
def rows_deg(gens,Dt):
    out=[]
    for g,D in gens:
        if D<=Dt: out+=E.ideal_rows(g,D,Dt)
    return out
