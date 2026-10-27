# strictness of one-modulus SL class survivors (inside eps-class mod n), first4 etc.
import sys, json, itertools, rlib, onemod
from math import gcd
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); F=json.load(open(sys.argv[3]))['F']; qs=json.loads(sys.argv[4])
for q0 in qs:
    pis,V,data=onemod.setup(P,K,F,q0); r=len(P)
    m,L=rlib.labels(pis,q0); n=m//4
    par=[x%2 for x in L]
    cls=set()
    for u in range(1,n+1):
        if gcd(u,n)!=1: continue
        for e in itertools.product((1,-1),repeat=r):
            cls.add(tuple((u*e[j]*L[j])%n for j in range(r)))
    U1=[u for u in range(1,m) if gcd(u,m)==1 and u%4==1]
    strict10=set(tuple((u*L[j])%m for j in range(r)) for u in U1)
    strictany=set()
    for e in itertools.product((1,-1),repeat=r):
        for t in itertools.product((0,1),repeat=r):
            base=[e[j]*L[j]+(m//2)*t[j] for j in range(r)]
            for u in U1: strictany.add(tuple((u*b)%m for b in base))
    surv=[]
    for c0 in cls:
        for k in itertools.product(range(4),repeat=r):
            c=tuple(c0[j]+n*k[j] for j in range(r))
            if any(c[j]%2!=par[j] for j in range(r)): continue
            ok=True
            for v in V:
                s=sum(a*b for a,b in zip(v,c))%m
                if s%n: continue
                D,N,g=data[v]; sp=s//n
                if q0*g[sp]>abs(D[sp]): ok=False; break
            if ok: surv.append(c)
    n10=sum(1 for c in surv if c in strict10); nany=sum(1 for c in surv if c in strictany)
    print(json.dumps({'q0':q0,'class_SL_survivors':len(surv),'strict_any':nany,'strict_10':n10}),flush=True)
