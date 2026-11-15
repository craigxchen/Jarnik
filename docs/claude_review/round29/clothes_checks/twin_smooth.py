# smooth-modulus twin: q = product of consecutive odd primes coprime to N, q >= lam*N^(1/4); max class size of z mod q
import sys, math, sympy, random
from math import gcd
from common import gpoints, arcmax
from twin import twin
def smooth_q(N, lam):
    q=1; Q=lam*N**0.25
    for l in sympy.primerange(3,10**6):
        if N%l: q*=l
        if q>=Q: return q
def classmax(P,q):
    d={}
    for z in P: d[(z[0]%q,z[1]%q)]=d.get((z[0]%q,z[1]%q),0)+1
    return max(d.values())
lam=float(sys.argv[1])
split=[p for p in sympy.primerange(5,2000) if p%4==1]
random.seed(1)
tests=[math.prod(split[:k]) for k in range(3,10)]
for _ in range(12):
    k=random.randint(4,8); ps=random.sample(split[:40],k)
    tests.append(math.prod(p**random.choice([1,1,1,2,3]) for p in ps))
for N in tests:
    P=gpoints(N); q=smooth_q(N,lam)
    (t,qq),r=twin(N,lam) if len(P)<=1100 else ((None,None),len(P))
    print(N,"r=",len(P),"smooth q=",q,"cls=",classmax(P,q)," best-any-odd-q:",t,qq," arc(C=2/lam):",arcmax(N,2/lam),flush=True)
