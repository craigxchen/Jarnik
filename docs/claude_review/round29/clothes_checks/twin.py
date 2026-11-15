# residue-class twin: T(N,lam)=max over q>=lam*N^(1/4), gcd(q,N)=1, of #{z: |z|^2=N, z==z0 mod q}
import sys, math, itertools, sympy
from math import gcd
from common import gpoints, arcmax
def twin(N, lam):
    P=gpoints(N); Q=lam*N**0.25
    best=(1,None)
    for z1 in P:
        cnt={}
        for z in P:
            if z==z1: continue
            c=gcd(z[0]-z1[0], z[1]-z1[1])
            for d in map(int,sympy.divisors(c)):
                if d>=Q and d%2==1 and gcd(d,N)==1: cnt[d]=cnt.get(d,0)+1
        for d,v in cnt.items():
            if v+1>best[0]: best=(v+1,d)
    return best, len(P)
if __name__=="__main__":
    lam=float(sys.argv[1]) if len(sys.argv)>1 else 1.0
    split=[p for p in sympy.primerange(5,400) if p%4==1]
    tests=[]
    for k in range(2,9): tests.append(math.prod(split[:k]))
    tests += [5**8, 5**4*13**4, 5**3*13**2*17**2, 5**2*13**2*17**2*29, 5*13*17*29*37*41*53*61*73]
    for N in tests:
        (t,q),r=twin(N,lam)
        C=2/lam
        print(N, "r=",r, "twin=",t,"q=",q, "N^1/4=%.1f"%N**0.25, "arc(C=%.2f)="%C, arcmax(N,C), flush=True)
