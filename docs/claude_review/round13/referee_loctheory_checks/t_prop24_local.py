# Is J(P)_e (P of Prop 2.4) generated mod grid by pullbacks from sub-projections with rho < 8 ?
import itertools, sys
from rlib import *
from echelon import Ech
b=5; j=2; D=7; m=8
G=list(itertools.product(range(-j,j+1),repeat=b))
Pset=[y for y in G if sum(y)%3==0]
cap=[2*j+1]*b
target=len(monos_eq(b,D,cap))-244   # = 11 (H_P(7)=244 exact)
print("target dim (J(P)/grid)_7 =",target)
def prim(v):
    from math import gcd
    g=0
    for x in v: g=gcd(g,abs(x))
    return g==1
dirs=[v for v in itertools.product(range(-1,2),repeat=b) if any(v) and prim(v)]
dirs=[v for v in dirs if [x for x in v if x][0]>0]
print("directions",len(dirs))
E=Ech()
def add_family(F):
    img=sorted(set(tuple(sum(a*c for a,c in zip(v,y)) for v in F) for y in Pset))
    k=len(F)
    # rho of image (exact)
    rho=None
    for e in range(0,40):
        if h_func_exact(img,e)==len(img): rho=e;break
    if rho>=m: return rho,0
    added=0
    for e in range(1,D+1):
        for g in J_space(img,e):
            pb=pullback(g,F,b)
            for mm in monos_eq(b,D-e,cap):
                r=reduce_cap({tuple(x+y for x,y in zip(kk,mm)):vv for kk,vv in pb.items()},cap)
                if r and E.add(r): added+=1
    return rho,added
for v in dirs:
    rho,a=add_family([v])
    if a: print("1-dir",v,"rho",rho,"added",a,"rank",E.rank())
print("after 1-dim families rank",E.rank()); sys.stdout.flush()
cnt=0
for v1,v2 in itertools.combinations(dirs,2):
    if E.rank()>=target: break
    rho,a=add_family([v1,v2]); cnt+=1
    if a: print("2-dir",v1,v2,"rho",rho,"added",a,"rank",E.rank()); sys.stdout.flush()
print("pairs tried",cnt,"final rank (mod p lower bound)",E.rank(),"target",target)
