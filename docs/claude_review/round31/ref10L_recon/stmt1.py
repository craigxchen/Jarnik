"""Brute-force test of B_int <= B_free <= 2 max(4,B_int) (needs a prime p > C^2),
and B_free <= max(4,B_int) at C=1, on small squarefree moduli, all delta vectors.
Also: for every max clique found, test that each p-side (p > C^2) of size>=5 with
EVERY admissible choice of magnitudes... (too many); we test a random choice."""
import sys, random
from math import gcd, isqrt
from itertools import product, combinations

def isqrt_floor_C(C2, q):  # floor(sqrt(C2*q))
    return isqrt(C2*q)

def max_clique(n, adj):
    best=[]
    def rec(R,P):
        nonlocal best
        if len(R)+len(P)<=len(best): return
        if not P:
            best=R[:]; return
        for v in list(P):
            if len(R)+len(P)<=len(best): return
            rec(R+[v],[w for w in P if w>v and adj[v][w]])
            P=[w for w in P if w!=v]
    rec([],list(range(n)))
    return best

def run(primes, C2, deltas, out):
    k=len(primes); q=1
    for p in primes: q*=p
    L=isqrt(C2*q)
    # allowed[pattern] = list of admissible magnitudes; pattern bit=1 means separated
    allowed={}
    for pat in range(1,1<<k):
        h=1
        for t,p in enumerate(primes):
            if not (pat>>t)&1: h*=p
        lst=[]
        for x in range(h, L+1, h):
            if gcd(x,q)!=h: continue
            if all(x%p in (deltas[t]%p,(-deltas[t])%p) for t,p in enumerate(primes) if (pat>>t)&1):
                lst.append(x)
        allowed[pat]=lst
    n=1<<k
    adj=[[bool(allowed.get(i^j)) if i!=j else False for j in range(n)] for i in range(n)]
    cl=max_clique(n,adj)
    Bfree=len(cl)
    Bint=0; bestY=None
    for eps in product((1,-1),repeat=k):
        Y=[y for y in range(0,L+1) if all(y%p in (0,(eps[t]*deltas[t])%p) for t,p in enumerate(primes))]
        if len(Y)>Bint: Bint=len(Y); bestY=Y
    big=any(p*1>C2 for p in primes)
    viol=False
    if Bint>Bfree: viol=True
    if big and Bfree>2*max(4,Bint): viol=True
    if C2==1 and Bfree>max(4,Bint): viol=True
    return Bfree,Bint,viol,cl,allowed

if __name__=="__main__":
    random.seed(1)
    sets=[[5,7,11],[7,11,13],[11,13,17],[13,17,19],[5,7,11,13],[7,11,13,17],[3,5,7,11],[5,7,11,13,17],[17,19,23,29],[11,13,17,19,23]]
    for primes in sets:
        for C2 in (1,2,4,9):
            q=1
            for p in primes:q*=p
            if C2*q>4*10**6: continue
            alld=list(product(*[range(1,(p-1)//2+1) for p in primes]))
            if len(alld)>400: alld=random.sample(alld,400)
            stats={}; nviol=0; worst=0
            for dl in alld:
                bf,bi,v,cl,al=run(primes,C2,dl,None)
                if v:
                    nviol+=1; print("VIOLATION",primes,C2,dl,bf,bi)
                stats[(bf,bi)]=stats.get((bf,bi),0)+1
                worst=max(worst,bf/max(4,bi))
            print(primes,"C^2=",C2,"cases",len(alld),"viol",nviol,"max Bfree/max(4,Bint)=%.3f"%worst, sorted(stats.items())[-6:], flush=True)
