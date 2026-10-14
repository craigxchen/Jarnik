# Prop 2.4 image: 3- and 4-dim projections (kernel spanned by small integer vectors): rho and contributions in degree 7
import itertools, sys
from rlib import *
from echelon import Ech
b=5; j=2; D=7; m=8
G=list(itertools.product(range(-j,j+1),repeat=b))
Pset=[y for y in G if sum(y)%3==0]
cap=[5]*b
E=Ech()
E.add(reduce_cap(pullback({(7,):1},[(1,1,1,1,1)],b),cap))
zs=[z for z in itertools.product((-1,0,1),repeat=b) if any(z) and [x for x in z if x][0]>0]
from collections import Counter
rhos=Counter()
def doK(Zrows):
    K,nul=flint.fmpz_mat(Zrows).nullspace()
    Fm=[[int(K[r,c]) for r in range(b)] for c in range(nul)]
    img=sorted(set(tuple(sum(a*x for a,x in zip(v,y)) for v in Fm) for y in Pset))
    rho=next(e for e in range(60) if h_func_exact(img,e)==len(img))
    rhos[(nul,rho)]+=1
    if rho>=m: return
    for e in range(1,D+1):
        for g in J_space(img,e):
            pb=pullback(g,Fm,b)
            for mm in monos_eq(b,D-e,cap):
                r=reduce_cap({tuple(x+y for x,y in zip(kk,mm)):vv for kk,vv in pb.items()},cap)
                if r: E.add(r)
for z in zs: doK([list(z)])
print("4-dim projections: (dim,rho) counts",sorted(rhos.items()),"rank",E.rank(),"/ 11"); sys.stdout.flush()
rhos.clear()
for z1,z2 in itertools.combinations(zs,2):
    if flint.fmpz_mat([list(z1),list(z2)]).rank()<2: continue
    doK([list(z1),list(z2)])
print("3-dim projections: (dim,rho) counts",sorted(rhos.items()),"rank",E.rank(),"/ 11")
