import math, numpy as np, json
from ref_seq import enum
for primes, T, shells in [([13,17,29,37,41,53,61,73],34,[20,25,28,31,34]), ([13,17,29,37,41,53,61,73,89,97],33,[23,27,31,33])]:
    logs=[math.log(p) for p in primes]
    Pi=np.prod([1+2/(math.sqrt(p)-1) for p in primes]); Lam=sum(p**-0.5 for p in primes)
    V=enum(logs,T); w=np.abs(V)@np.array(logs)
    tot=(Pi-1)/2  # sum over v != 0 modulo sign of e^{-w/2}
    s=lambda a,b: float(np.exp(-w[(w>a)&(w<=b)]/2).sum())
    sh={f"({shells[i]},{shells[i+1]}]": round(2.2*s(shells[i],shells[i+1]),4) for i in range(len(shells)-1)}
    tail=2.2*(tot-float(np.exp(-w/2).sum()))
    print(json.dumps(dict(primes=primes,Lambda=round(Lam,4),Pi=round(Pi,3),nV=len(V),shells=sh,tail_beyond_T=round(tail,4),total_expected=round(2.2*tot,2))))
