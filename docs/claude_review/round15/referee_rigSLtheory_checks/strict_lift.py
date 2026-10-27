import sys, json, itertools, random, rlib, onemod
from math import gcd
P=[5,13,17,29]; K=int(sys.argv[1]); F=json.load(open(sys.argv[2]))['F']; q0=int(sys.argv[3])
pis,V,data=onemod.setup(P,K,F,q0); r=4
m,L=rlib.labels(pis,q0); n=m//4; par=[x%2 for x in L]
U1=[u for u in range(1,m) if gcd(u,m)==1 and u%4==1]
strictany=set()
for e in itertools.product((1,-1),repeat=r):
    for t in itertools.product((0,1),repeat=r):
        base=[e[j]*L[j]+(m//2)*t[j] for j in range(r)]
        for u in U1: strictany.add(tuple((u*b)%m for b in base))
cls=set()
for u in range(1,n+1):
    if gcd(u,n)!=1: continue
    for e in itertools.product((1,-1),repeat=r): cls.add(tuple((u*e[j]*L[j])%n for j in range(r)))
rng=random.Random(3); found=0
for c0 in sorted(cls):
    for k in itertools.product(range(4),repeat=r):
        c=tuple(c0[j]+n*k[j] for j in range(r))
        if any(c[j]%2!=par[j] for j in range(r)) or c in strictany: continue
        ok=all(not (sum(a*b for a,b in zip(v,c))%m%n==0 and q0*data[v][2][(sum(a*b for a,b in zip(v,c))%m)//n]>abs(data[v][0][(sum(a*b for a,b in zip(v,c))%m)//n])) for v in V)
        if not ok: continue
        # level-2 lift without level-2 collisions
        m2=q0*m; lift=None
        for trial in range(20000):
            t=[rng.randrange(q0) for _ in range(r)]; cc=[c[j]+m*t[j] for j in range(r)]
            if all(sum(a*b for a,b in zip(v,cc))%(m2//4)!=0 for v in V): lift=cc; break
        print('non-strict SL label',c,'lift',lift); found+=1
        if found>=3: sys.exit()
