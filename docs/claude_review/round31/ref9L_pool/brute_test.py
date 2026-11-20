import random, itertools, sys
from math import isqrt, prod
sys.path.insert(0,'.')
from twin_sim import max_clique, admissible_triples
def brute(primes,L):
    q=prod(primes); best=[0]
    def adm(S):
        return all(len({y%p for y in S})<=2 for p in primes)
    def rec(cur,start):
        nonlocal best
        if len(cur)>len(best): best=list(cur)
        for c in range(start,L+1):
            if adm(cur+[c]): rec(cur+[c],c+1)
    rec([0],1); return best
rng=random.Random(1)
from sympy import primerange
for t in range(40):
    pr=sorted(rng.sample(list(primerange(5,40)),3)); q=prod(pr)
    for C in (1,2,3):
        L=isqrt(C*C*q)
        b1=brute(pr,L); b2,T=max_clique(pr,L,q)
        Tb={(a,b) for a in range(1,L+1) for b in range(a+1,L+1) if (a*b*(b-a))%q==0}
        assert Tb==T,(pr,C)
        assert len(b1)==len(b2),(pr,C,b1,b2)
print('ok')
