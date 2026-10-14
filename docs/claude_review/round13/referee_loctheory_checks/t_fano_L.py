import itertools, sys, time, pickle
from rlib import *
from echelon import Ech
Y,h,dJ6,dJ7=pickle.load(open('fano_Y.pkl','rb'))
j=2; cap=[5]*7
avec=[a for a in itertools.product((0,1),repeat=3) if any(a)]
zs=list(itertools.product((0,1),repeat=3))
dot=lambda a,z: sum(x*y for x,y in zip(a,z))%2
idx={a:i for i,a in enumerate(avec)}
add=lambda a,b: tuple((x+y)%2 for x,y in zip(a,b))
lines=set(); quads=set()
for a,b in itertools.combinations(avec,2):
    lines.add(frozenset([a,b,add(a,b)]))
for a,b,c in itertools.combinations(avec,3):
    d=add(add(a,b),c)
    if any(d) and d not in (a,b,c) and add(a,b)!=c: quads.add(frozenset([a,b,c,d]))
print(len(lines),"lines",len(quads),"quads")
gens={}   # degree -> list of homogeneous dicts in 7 vars
def addgen(e,g): gens.setdefault(e,[]).append(g)
for Ls in list(lines)+list(quads):
    ids=sorted(idx[a] for a in Ls); c=len(ids)
    img=sorted(set(tuple(y[i] for i in ids) for y in Y))
    for e in range(1,8):
        for g in J_space(img,e,[5]*c):
            addgen(e,{tuple(ex[ids.index(q)] if q in ids else 0 for q in range(7)):v for ex,v in g.items()})
print("line/quad local gens by degree",{e:len(v) for e,v in gens.items()})
# width powers of class-sum combinations: l_m(y) = sum_z m_z c_z(y), c_z = -(1/4) sum_a y_a chi_a(z); use 4*l_m (integer form)
seen=set(); nw=0
for m in itertools.product((-1,0,1),repeat=8):
    coef=[-sum(mz*(-1)**dot(a,z) for mz,z in zip(m,zs)) for a in avec]
    if not any(coef): continue
    from math import gcd
    g=0
    for x in coef: g=gcd(g,abs(x))
    coef=tuple(x//g for x in coef)
    if [x for x in coef if x][0]<0: coef=tuple(-x for x in coef)
    if coef in seen: continue
    seen.add(coef)
    vals=set(sum(cc*yy for cc,yy in zip(coef,y)) for y in Y)
    if len(vals)<=7:
        f=pullback({(len(vals),):1},[coef],7); addgen(len(vals),f); nw+=1
print("distinct primitive class-combination directions",len(seen),"with <=7 values:",nw)
for D in (6,7):
    rows=ideal_degree_rows(gens,7,D,cap)
    cols=monos_eq(7,D,cap)
    rp=rank_rows(rows,cols,modp=True)
    print("D",D,"rows",len(rows),"rank mod p",rp, "target",dJ6 if D==6 else dJ7); sys.stdout.flush()
    if D==7:
        T0=time.time(); rq=rank_rows(rows,cols,modp=False); print("exact Q rank",rq,time.time()-T0)
pickle.dump(gens,open('fano_gens.pkl','wb'))
