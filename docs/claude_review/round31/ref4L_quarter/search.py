import sys, random, gmpy2
from itertools import combinations
from math import prod
x=(-27,-270,135,120,40)
s=[sum(prod(c) for c in combinations(x,k)) for k in range(6)]
E=list(combinations(range(5),2))
b=[x[i]*x[j]*(s[1]-x[i]-x[j]) for i,j in E]
q=int(sys.argv[2]); D=int(sys.argv[1]); seed=int(sys.argv[3])
random.seed(seed)
small=[p for p in range(7,2000) if gmpy2.is_prime(p)]
start=random.randrange(10**D,2*10**D); start-= start%60; start+=1  # P = 1 mod 60 -> coprime to 30, 1 mod 4
# allowed residues
bad={p:set((-q*bb)%p for bb in b) for p in small}
P=start; tries=0
while True:
    P+=60
    if any(P%p in bad[p] for p in small): continue
    tries+=1
    if all(gmpy2.is_prime(P+q*bb, 30) for bb in b):
        print(P, tries); break
