# AG(3,2) image at j=1, degrees 4..6: do 7-direction projections (kernel z in {-1,0,1}^8) close the deficit over
# grid + circuits + width powers?  Exact generators, mod-p rank vs exact target (sandwich).
import itertools, sys, flint, pickle
from rlib import *
from echelon import Ech
j=1; b=8; t=4
V=[a for a in itertools.product((0,1),repeat=t) if any(a)]
F=[a for a in V if sum(a)%2==1]
codes=[s for s in itertools.product((0,1),repeat=b) if any(s) and all(sum(F[i][q] for i in range(b) if s[i])%2==0 for q in range(t))]
P=[y for y in itertools.product(range(-j,j+1),repeat=b) if sum(y)%4==0 and all(sum(y[i] for i in range(b) if s[i])%2==0 for s in codes)]
cap=[3]*b
circ=[s for s in codes if not any(u!=s and all(ui<=si for ui,si in zip(u,s)) for u in codes)]
target={4:182,5:504,6:784}
gens={}
for s in circ:
    idx=[q for q in range(b) if s[q]]
    Ps=[y for y in itertools.product(range(-j,j+1),repeat=4) if sum(y)%2==0]
    for e in range(1,7):
        for g in J_space(Ps,e,[3]*4):
            gens.setdefault(e,[]).append({tuple(ex[idx.index(q)] if q in idx else 0 for q in range(b)):v for ex,v in g.items()})
for v in itertools.product((-1,0,1),repeat=b):
    if not any(v) or [x for x in v if x][0]<0: continue
    vals=set(sum(a*x for a,x in zip(v,y)) for y in P)
    if len(vals)<=6: gens.setdefault(len(vals),[]).append(pullback({(len(vals),):1},[v],b))
for D in (4,5,6):
    E=Ech()
    for r in ideal_degree_rows(gens,b,D,cap): E.add(r)
    start=E.rank(); used=0
    for z in itertools.product((-1,0,1),repeat=b):
        if not any(z) or [x for x in z if x][0]<0: continue
        if E.rank()>=target[D]: break
        K,nul=flint.fmpz_mat([list(z)]).nullspace()
        Fm=[[int(K[r,c]) for r in range(b)] for c in range(nul)]
        img=sorted(set(tuple(sum(a*x for a,x in zip(v,y)) for v in Fm) for y in P))
        before=E.rank()
        ng={}
        for e in range(1,D+1):
            sp=J_space(img,e)
            if sp: ng[e]=[pullback(g,Fm,b) for g in sp]
        for r in ideal_degree_rows(ng,b,D,cap): E.add(r)
        if E.rank()>before: used+=1
    print("D",D,"target",target[D],"circuits+width",start,"+ 7-dir projections (z in {-1,0,1}^8)",E.rank(),"families used",used); sys.stdout.flush()
