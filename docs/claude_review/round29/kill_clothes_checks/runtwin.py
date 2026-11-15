import sys, math, random, sympy, kc
lam=float(sys.argv[1])
s1=[p for p in sympy.primerange(5,400) if p%4==1]
A=[math.prod(s1[:k]) for k in range(2,9)]+[5**8,5**4*13**4,5**3*13**2*17**2,5**2*13**2*17**2*29,5*13*17*29*37*41*53*61*73]
s2=[p for p in sympy.primerange(5,2000) if p%4==1]; random.seed(1)
B=[math.prod(s2[:k]) for k in range(3,10)]
for _ in range(12):
    k=random.randint(4,8); ps=random.sample(s2[:40],k)
    B.append(math.prod(p**random.choice([1,1,1,2,3]) for p in ps))
L=sorted(set(A+B))
for N in L:
    P=kc.points(N)
    t=kc.twinmax(N,lam,P) if len(P)<=1100 else (None,)
    print(N,"inA" if N in A else "", "inB" if N in B else "", "r=",len(P),"twin=",t[:2],"arc(C=2/lam)=",kc.arcmax(N,2/lam),flush=True)
