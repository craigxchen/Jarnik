import numpy as np, itertools, pickle, sys
from flint import nmod_mat, fmpz_mat
from paley_base import build, COLS
from reng import REng, PR
d=11; Dm=7
A=build(2)
shift=[(j+1)%11 for j in range(11)]   # (sigma f)(x_j) = x_{j+1}
# point orbits
seen=set(); reps=[]
for k in map(tuple,A):
    if k in seen: continue
    orb=set(); cur=k
    while cur not in orb:
        orb.add(cur); cur=cur[1:]+cur[:1]
    seen|=orb; reps.append(k)
print("point orbits",len(reps))
E=REng(d,Dm,shift)
K=np.array(reps,dtype=np.int64)
Zr=(K==0).astype(np.int64); Nr=(K==-1).astype(np.int64)
cols=[]; degs=[]
for D in range(Dm+1):
    M=np.array(E.ms[D],dtype=np.int64)
    zero=(Zr@(M>0).T.astype(np.int64))>0
    par=(Nr@(M%2==1).T.astype(np.int64))%2
    val=np.where(zero,0,1-2*par)  # points x monomials
    oid,cnt=E.orb[D]
    S=np.zeros((len(reps),cnt),dtype=np.int64)
    for o in range(cnt): pass
    np.add.at(S.T, oid, val.T)
    cols.append(S); degs+= [D]*cnt
    print("deg",D,"monomials",len(E.ms[D]),"orbits",cnt)
Ecat=np.concatenate(cols,axis=1)
print("E inv shape",Ecat.shape)
np.save("paley_Einv.npy",Ecat); pickle.dump(degs,open("paley_Einv_degs.pkl","wb"))
# mod p cumulative ranks -> upper bounds on dim P0 J_D
cum=[]
c0=0
prev=0
for D in range(Dm+1):
    c0+=cols[D].shape[1]
    sub=Ecat[:, :c0]
    M=nmod_mat(sub.shape[0],sub.shape[1],[int(x)%PR for x in sub.flatten()],PR)
    r=M.rank()
    HD=r-prev; prev=r
    print("D",D,"inv H",HD,"inv dim J_D upper bound (mod p)",cols[D].shape[1]-HD)
