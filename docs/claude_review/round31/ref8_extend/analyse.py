# For every recorded genuine window: normalise, exponent affine deficiency def0 = M-1-affrank, and the shortest
# prefix S_k = {3,7,11,...}(first k inert primes) whose residue indicators restore affine rank M-1.
# Exact: rank mod P61 is a lower bound for rank over Q; when it is not full we recompute over Q exactly.
import pickle, sys, math, collections
from gaussx import *
from math import log
P61=(1<<61)-1
from sympy import primerange
INERT=[l for l in primerange(3,5000) if l%4==3]
def rank_mod(rows,p=P61):
    A=[[x%p for x in r] for r in rows]; m=len(A); n=len(A[0]) if A else 0; r=0
    for c in range(n):
        piv=next((i for i in range(r,m) if A[i][c]),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]; inv=pow(A[r][c],p-2,p)
        for i in range(m):
            if i!=r and A[i][c]:
                f=A[i][c]*inv%p; A[i]=[(x-f*y)%p for x,y in zip(A[i],A[r])]
        r+=1
        if r==m: break
    return r
def affrank(rows):
    if len(rows)<=1: return 0
    D=[[x-y for x,y in zip(r,rows[0])] for r in rows[1:]]
    if not D[0]: return 0
    rm=rank_mod(D)
    if rm==len(D): return rm
    return rank_Q(D)
def sharpK(C):
    t=log(C)-0.5*log(2); s=0; k=0
    while s<=t: s+=log(INERT[k])/(INERT[k]+1); k+=1
    return k
def analyse(ps,es,N,res):
    out=[]
    for M,C,zs,levs in res:
        cols=list(zip(*levs))
        mins=[min(c) for c in cols]; maxs=[max(c) for c in cols]
        keep=[k for k in range(len(ps)) if maxs[k]>mins[k]]
        rows=[[levs[i][k]-mins[k] for k in keep] for i in range(M)]
        om=len(keep)
        if om==0: continue
        ar=affrank(rows); d0=M-1-ar
        need=0
        if d0>0:
            # normalised marks z/g
            g=(1,0)
            for k,p in enumerate(ps):
                g=gmul(g,gmul(gpow(PI[p],mins[k]),gpow(gconj(PI[p]),es[k]-maxs[k])))
            zz=[gdivexact(z,g) for z in zs]
            ext=[list(r) for r in rows]; k=0
            while affrank(ext)<M-1:
                l=INERT[k]; cls=sorted(set((z[0]%l,z[1]%l) for z in zz))
                for i in range(M): ext[i]+= [1 if (zz[i][0]%l,zz[i][1]%l)==c else 0 for c in cls]
                k+=1
            need=k
        Wp=sum((maxs[k]-mins[k])*log(ps[k]) for k in keep)
        out.append((M,C,om,d0,need,Wp,N))
    return out
from enum_clusters import PI
if __name__=="__main__":
    allr=[]
    for f in sys.argv[1:]:
        for ps,es,N,res in pickle.load(open(f,"rb")):
            allr+=analyse(ps,es,N,res)
    pickle.dump(allr,open("analysed.pkl","wb"))
    print("windows analysed",len(allr))
    # 1. VIII.1 sharpened: need <= sharpK(C) always (theorem).  Also C<=sqrt2 => d0=0.
    viol=[r for r in allr if r[4]>sharpK(max(r[1],1e-9))]
    print("violations of sharpened VIII.1 (need > |S_C|):",len(viol))
    print("windows with C<=sqrt2 and d0>0:",sum(1 for r in allr if r[1]<=2**.5 and r[3]>0))
    # 2. smallest C with d0>=k
    for k in range(1,8):
        c=[r for r in allr if r[3]>=k]
        if c:
            b=min(c,key=lambda r:r[1]); print(f"def0>={k}: min C={b[1]:.4f} (M={b[0]}, omega'={b[2]}, N={b[6]})  count={len(c)}")
    # 3. by C-bins: max d0, max need, and max (M-1-omega')
    bins=[1,1.5,2,3,4,5,6,8,10,12,16]
    for lo,hi in zip(bins,bins[1:]):
        s=[r for r in allr if lo<r[1]<=hi]
        if s: print(f"C in ({lo},{hi}]: n={len(s)} maxM={max(r[0] for r in s)} max def0={max(r[3] for r in s)} max need |S|={max(r[4] for r in s)} (bound {sharpK(hi)}) max(M-1-omega')={max(r[0]-1-r[2] for r in s)}")
