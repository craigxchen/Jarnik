import sys, random, gmpy2, numpy as np
from itertools import combinations
from math import prod
x=(-27,-270,135,120,40)
s=[sum(prod(c) for c in combinations(x,k)) for k in range(6)]
E=list(combinations(range(5),2))
b=[x[i]*x[j]*(s[1]-x[i]-x[j]) for i,j in E]
D=int(sys.argv[1]); q=int(sys.argv[2]); random.seed(int(sys.argv[3]))
small=[p for p in range(7,int(sys.argv[4])) if gmpy2.is_prime(p)]
W=1<<22
start=random.randrange(10**D,2*10**D); start-=start%60; start+=1
tested=0
while True:
    ok=np.ones(W,dtype=bool)
    for p in small:
        inv60=pow(60,-1,p)
        for bb in b:
            # P=start+60k, need start+60k+q*bb != 0 mod p
            r=(-(start+q*bb))*inv60 % p
            ok[r::p]=False
    for kk in np.nonzero(ok)[0]:
        P=start+60*int(kk); tested+=1
        if all(gmpy2.is_prime(P+q*bb,30) for bb in b):
            print(P,q,tested,flush=True); sys.exit()
    start+=60*W
