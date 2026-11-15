import math, random, sympy, kc
s2=[p for p in sympy.primerange(5,2000) if p%4==1]; random.seed(1)
B=[math.prod(s2[:k]) for k in range(3,10)]
for _ in range(12):
    k=random.randint(4,8); ps=random.sample(s2[:40],k)
    B.append(math.prod(p**random.choice([1,1,1,2,3]) for p in ps))
for lam in [1.0,0.5]:
    out=[]
    for N in B:
        P=kc.points(N); Q=lam*N**0.25; q=1
        for l in sympy.primerange(3,10**6):
            if N%l: q*=l
            if q>=Q: break
        out.append(max(len(c) for c in kc.classes(P,q).values()))
    print("lam",lam,"smooth-q max class per norm:",out,"max",max(out))
