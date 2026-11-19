"""Independent construction of the V.2 full-isotropy model on my own base parent
(four-block star of V.3 at my own parameter P=4445280047, D=6350400).
Checks: CRT SL2 lift, shear, all 20 endpoint incidences, signed cofactors preserved,
width ~ m^-2, genuine height H via Gaussian gcd, Liouville/CC consistency,
cycle flip -> all 25 balanced incidences."""
import math, sys
from itertools import combinations
from math import gcd, prod
from sympy import isprime, sqrt_mod
from sympy.ntheory.modular import crt
from gauss import *
exec(open('base_search.py').read().split("D=6350400")[0])
D=6350400; P=4445280047
cols,blocks,cof=model(D,P)
factors=[[140057,1282121,110061697],[97,203737018017362497],[541817,30557,1193544061],[61,757,427943780183617],
 [57345697,347376976897],[633877,6704653,4652689],[1061,18673684486469069],[29,433,1574294728196837],
 [17,1164243030535397777],[229,557,6959881,22261913]]
for f,J in zip(factors,blocks):
    assert prod(f)==J and all(isprime(p) and p%4==1 for p in f)
allp=[p for f in factors for p in f]; assert len(set(allp))==len(allp)
N=prod(blocks); print("parent digits",len(str(N)),"primes",len(allp))
edges=list(combinations(range(5),2))
for (i,j),J,t in zip(edges,blocks,cof):
    assert det(cols[i],cols[j])==t*J and gcd(abs(t),N)==1
assert all(gcd(*c)==1 for c in cols)
iso=[(i,e) for e,(a,b) in enumerate(edges) for i in (a,b) if norm(cols[i])%blocks[e]==0]
print("base model: isotropic endpoint incidences",len(iso),"of 20; blocks:",sorted(set(e for _,e in iso)))
def rel(x): return math.log10(x) if x>0 else -999
def width(V):
    # max doubled chord = max 2|det|/(|V_i||V_j|), in log10
    return max(math.log10(2*abs(det(V[i],V[j])))-0.5*(math.log10(norm(V[i]))+math.log10(norm(V[j]))) for i,j in edges)
print("base width log10 = %.2f ; N^-3/10 log10 = %.2f"%(width(cols), -0.3*math.log10(N)))
# marked line at p: column of endpoint i of the block (both endpoints on the same line mod p)
lres=[];ures=[];mods=[]
for e,(a,b) in enumerate(edges):
    for p in factors[e]:
        x,y=cols[a][0]%p, cols[a][1]%p
        xb,yb=cols[b][0]%p, cols[b][1]%p
        assert (x*yb-y*xb)%p==0
        l=1 if y==0 else 0
        y2=(y+l*x)%p; assert y2
        ip=sqrt_mod(-1,p)
        u=((ip*y2-x)*pow(y2,-1,p))%p
        lres.append(l); ures.append(u); mods.append(p)
Lc=int(crt(mods,lres)[0]); Uc=int(crt(mods,ures)[0])
def apply(V,l,u):
    out=[]
    for x,y in V:
        y=y+l*x; x=x+u*y; out.append((x,y))
    return out
def check_endpoint(V):
    assert all(gcd(*v)==1 for v in V)
    for e,(a,b) in enumerate(edges):
        J=blocks[e]
        assert norm(V[a])%J==0 and norm(V[b])%J==0 and dot(V[a],V[b])%J==0
        assert det(V[a],V[b])==cof[e]*J
    return True
res={}
for m in [0,1,10**6,N,N**3]:
    V=apply(cols,Lc,Uc+m*N)
    check_endpoint(V)
    zs,H=genuine(V) if m in (0,1,10**6,N) else (None,None)
    w=width(V)
    msg="m=%s: all 20 endpoint incidences isotropic, cofactors fixed; log10 width=%.1f"%("N^3" if m==N**3 else ("N" if m==N else m),w)
    if H:
        msg+="; log10 H=%.1f, log10 N=%.1f, log10(H^-3/10)=%.1f, log10 c_i ~ %s"%(math.log10(H),math.log10(N),-0.3*math.log10(H),
             [round(math.log10(norm(v))-sum(math.log10(blocks[e]) for e,(a,b) in enumerate(edges) if i in (a,b)),1) for i,v in enumerate(V)])
        # Cilleruelo-Cordoba type consistency: width must not beat H^{-3/10} by more than constants
        assert w >= -0.3*math.log10(H)-3
    print(msg); res[m]=V
# cycle flip on the m=N model
V=list(res[N]); sides=[set(e) for e in edges]
nn=lambda V,i: prod(blocks[e] for e in range(10) if i in sides[e])
ratio0=[ (norm(V[i]), nn(V,i)) for i in range(5)]
for e,(a,b) in enumerate(edges):
    if not ((b-a)==1 or (a,b)==(0,4)): continue
    J=blocks[e]
    g=ggcd((J,0),V[a]); assert norm(g)==J
    assert divides(g,V[b])
    gb=conj(g); newV=[]
    for k,v in enumerate(V):
        w_=mul(v,gb)
        if k in sides[e]:
            assert w_[0]%J==0 and w_[1]%J==0; w_=(w_[0]//J,w_[1]//J)
        newV.append(w_)
    V=newV; sides[e]=set(range(5))-sides[e]
assert all(sum(i in S for S in sides)==5 for i in range(5))
assert all(gcd(*v)==1 for v in V)
for S,J in zip(sides,blocks):
    for i in S: assert norm(V[i])%J==0
for e,(i,j) in enumerate(edges):
    g=prod(blocks[f] for f in range(10) if {i,j}<=sides[f])
    d=det(V[i],V[j]); assert d%g==0 and d//g==cof[e] and dot(V[i],V[j])%g==0
for i in range(5):
    assert norm(V[i])*ratio0[i][1]==ratio0[i][0]*nn(V,i)
zs2,H2=genuine(V); zs1,H1=genuine(res[N])
# same doubled ratio point: zeta' proportional to zeta
q=[mul(zs2[i],conj(zs1[i])) for i in range(5)]
assert all(q[i][0]*q[0][1]==q[i][1]*q[0][0] for i in range(5)) and H1==H2
print("cycle flip: all 25 balanced incidences, cofactors and |V_i|^2/n_i preserved, same genuine tuple (H equal)")
