# EVIDENCE ONLY (J spaces of images mod p). Do 6-dimensional projections (families of 6 directions, kernel z in {-1,0,1}^7) of Y_3(2) close the degree-7 gap?
import itertools, sys, time, pickle
from math import gcd
from rlib import *
from echelon import Ech
Y,h,dJ6,dJ7=pickle.load(open('fano_Y.pkl','rb'))
gens=pickle.load(open('fano_gens.pkl','rb'))
cap=[5]*7; D=7
avec=[a for a in itertools.product((0,1),repeat=3) if any(a)]
# GL(3,2) as coordinate permutations
def mat_apply(Mx,a): return tuple(sum(Mx[r][c]*a[c] for c in range(3))%2 for r in range(3))
perms=[]
for Mx in itertools.product(itertools.product((0,1),repeat=3),repeat=3):
    imgs=[mat_apply(Mx,a) for a in avec]
    if len(set(imgs))==7 and all(any(x) for x in imgs):
        perms.append([avec.index(x) for x in imgs])
print("group size",len(perms))
def permpoly(f,p):   # variable i -> variable p[i]
    out={}
    for ex,v in f.items():
        ne=[0]*7
        for i,a in enumerate(ex): ne[p[i]]=a
        out[tuple(ne)]=v
    return out
def J_space_p(X,e,p=2147483647):
    ms=monos_le(len(X[0]),e)
    M=flint.nmod_mat([[evalmono(q,m)%p for m in ms] for q in X],p)
    N,nul=M.nullspace()
    top=[i for i,m in enumerate(ms) if sum(m)==e]
    if nul==0: return []
    R=flint.nmod_mat([[int(N[i,c]) for i in top] for c in range(nul)],p)
    Rr,rk=R.rref()
    out=[]
    for r in range(rk):
        out.append({ms[top[c]]:int(Rr[r,c]) for c in range(len(top)) if int(Rr[r,c])})
    return out
E=Ech()
for D_ in [D]:
    for r in ideal_degree_rows(gens,7,D,cap): E.add(r)
print("start rank",E.rank(),"target",dJ7); sys.stdout.flush()
# orbit representatives of z under GL(3,2) x sign
zs=[z for z in itertools.product((-1,0,1),repeat=7) if any(z)]
seen=set(); reps=[]
for z in zs:
    if z in seen: continue
    orb=set()
    for p in perms:
        w=[0]*7
        for i in range(7): w[p[i]]=z[i]
        orb.add(tuple(w)); orb.add(tuple(-x for x in w))
    seen|=orb; reps.append(z)
print("orbit reps",len(reps)); sys.stdout.flush()
reps.sort(key=lambda z: sum(abs(x) for x in z))
for z in reps:
    T0=time.time()
    K,nul=flint.fmpz_mat([list(z)]).nullspace()
    F=[[int(K[r,c]) for r in range(7)] for c in range(nul)]
    img=sorted(set(tuple(sum(a*b for a,b in zip(v,y)) for v in F) for y in Y))
    newg={}
    for e in range(1,D+1):
        sp=J_space_p(img,e)
        if sp: newg[e]=[pullback(g,F,7) for g in sp]
    before=E.rank()
    for p in perms:
        pg={e:[permpoly(g,p) for g in gl] for e,gl in newg.items()}
        for r in ideal_degree_rows(pg,7,D,cap): E.add(r)
    print("z",z,"|img|",len(img),"rank",before,"->",E.rank(),"/",dJ7,"%.0fs"%(time.time()-T0)); sys.stdout.flush()
    if E.rank()>=dJ7: print("CLOSED"); break
