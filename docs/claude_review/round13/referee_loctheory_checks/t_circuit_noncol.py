# Are the degree-jc generators of a circuit image P_s = {y in [-j,j]^c : sum y even} generated (mod grid)
# by pullbacks along families of < c directions (non-column directions allowed)?
import itertools, sys
from math import gcd
from rlib import *
from echelon import Ech
c=int(sys.argv[1]); j=int(sys.argv[2]); R=int(sys.argv[3]); kmax=int(sys.argv[4])
D=j*c
G=list(itertools.product(range(-j,j+1),repeat=c))
Pset=[y for y in G if sum(y)%2==0]
cap=[2*j+1]*c
HG=len(monos_eq(c,D,cap)); target=(HG-(-1)**(D))//2
print("c",c,"j",j,"degree",D,"target dim (J(P_s)/grid)_D =",target,"(Thm 3.3 count)")
def prim(v):
    g=0
    for x in v: g=gcd(g,abs(x))
    return g==1
dirs=[v for v in itertools.product(range(-R,R+1),repeat=c) if any(v) and prim(v)]
dirs=[v for v in dirs if [x for x in v if x][0]>0]
E=Ech(); maxrho=0; used=[]
def add_family(F):
    global maxrho
    img=sorted(set(tuple(sum(a*x for a,x in zip(v,y)) for v in F) for y in Pset))
    added=0
    for e in range(1,D+1):
        for g in J_space(img,e):
            pb=pullback(g,F,c)
            for mm in monos_eq(c,D-e,cap):
                r=reduce_cap({tuple(x+y for x,y in zip(kk,mm)):vv for kk,vv in pb.items()},cap)
                if r and E.add(r): added+=1
    if added:
        rho=next(e for e in range(200) if h_func_exact(img,e)==len(img))
        maxrho=max(maxrho,rho); used.append((F,rho,added))
    return added
import random
for k in range(1,kmax+1):
    fams=list(itertools.combinations(dirs,k))
    random.seed(1); random.shuffle(fams)
    for F in fams:
        if E.rank()>=target: break
        M=flint.fmpz_mat([list(v) for v in F])
        if M.rank()<k: continue
        add_family(list(F))
    print("after families of <=",k,"directions (entries in [-%d,%d]): rank mod p ="%(R,R),E.rank(),"/",target); sys.stdout.flush()
    if E.rank()>=target: break
print("families used",len(used),"max rho among used",maxrho)
for u in used[:12]: print(u)
