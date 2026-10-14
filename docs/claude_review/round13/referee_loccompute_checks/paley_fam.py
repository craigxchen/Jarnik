# Incremental mod-p rank of P0(L_7) for Paley-12 L=2 with arbitrary integer-direction families.
import numpy as np, itertools, sys, pickle, time
from flint import nmod_mat
from paley_base import build, COLS
from reng import REng, PR
from rlib import J_exact, monos
d=11; DT=7
A=build(2)
shift=[(j+1)%11 for j in range(11)]
E=REng(d,DT,shift)
class Acc:
    def __init__(self):
        self.B=np.zeros((0,E.orb[DT][1]),dtype=np.int64); self.r=0
    def add(self, rows):
        if not rows: return 0
        M=np.vstack([self.B]+[np.array(rows,dtype=np.int64)])
        nm=nmod_mat(M.shape[0],M.shape[1],[int(x) for x in M.flatten()],PR)
        R,rk=nm.rref()
        newr=rk-self.r
        if newr>0:
            self.B=np.array([[int(R[i,j]) for j in range(M.shape[1])] for i in range(rk)],dtype=np.int64)
            self.r=rk
        return newr
def gen_rows(g, D):
    return [E.P0(v,DT) for v in E.ideal_rows(g,D,DT)]
def width_rows(v):
    vals=np.unique(A@v); s=len(vals)
    if s>DT: return []
    g=E.pullback([1],[(s,)],[v])
    return gen_rows(g,s)
def family_rows(U, Dmax=DT, skip_grid=True):
    U=[np.array(u,dtype=np.int64) for u in U]
    img=np.unique(np.stack([A@u for u in U],axis=1),axis=0)
    sizes=[len(np.unique(img[:,i])) for i in range(len(U))]
    if skip_grid and np.prod(sizes)==len(img): return None, len(img)
    J=J_exact([tuple(int(x) for x in p) for p in img], Dmax)
    rows=[]
    for D in range(Dmax+1):
        tm=monos(len(U),D)
        for vec in J[D]:
            g=E.pullback(vec,tm,U)
            if g.any(): rows+=gen_rows(g,D)
    return rows, len(img)
