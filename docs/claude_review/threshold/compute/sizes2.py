import sys
from math import factorial
from functools import lru_cache
sys.setrecursionlimit(10000)
def partitions(n, maxpart=None, maxlen=None):
    if maxpart is None: maxpart=n
    if n==0: yield (); return
    if maxlen==0: return
    for k in range(min(n,maxpart),0,-1):
        for rest in partitions(n-k,k,None if maxlen is None else maxlen-1):
            yield (k,)+rest
def orbits(M,p,n):
    out=[]
    Ps = list(partitions(p,maxlen=M)) if p>0 else [()]
    for P in Ps:
        Ns = list(partitions(n,maxlen=M-len(P))) if n>0 else [()]
        for Nn in Ns:
            z=M-len(P)-len(Nn)
            v=tuple(sorted(list(P)+[-x for x in Nn]+[0]*z,reverse=True))
            out.append(v)
    return sorted(set(out))
def mu_of(v):
    from collections import Counter
    return tuple(sorted(Counter(v).values(),reverse=True))
@lru_cache(None)
def kostka(lam, mu):
    # number of SSYT shape lam content mu (mu as tuple of multiplicities, order irrelevant)
    # recursive: remove horizontal strip of size mu[-1]
    if sum(lam)!=sum(mu): return 0
    if not mu: return 1
    m=mu[-1]; rest=mu[:-1]
    tot=0
    # remove horizontal strip of size m from lam: new shape nu with lam/nu horizontal strip
    L=len(lam)
    def rec(i, rem, nu):
        nonlocal tot
        if i==L:
            if rem==0:
                nu2=tuple(x for x in nu if x>0)
                tot_add=kostka(nu2, rest)
                return tot_add
            return 0
        s=0
        lo = lam[i+1] if i+1<L else 0
        for r in range(0, min(rem, lam[i]-lo)+1):
            s+=rec(i+1, rem-r, nu+[lam[i]-r])
        return s
    return rec(0,m,[])
def flam(lam):
    n=sum(lam); conj=[sum(1 for x in lam if x>j) for j in range(lam[0])]
    h=1
    for i,r in enumerate(lam):
        for j in range(r): h*= (r-j-1)+(conj[j]-i-1)+1
    return factorial(n)//h
if __name__=='__main__':
    M=int(sys.argv[1]); p=int(sys.argv[2]); n=int(sys.argv[3])
    O=orbits(M,p,n)
    mus=[mu_of(v) for v in O]
    tot=0
    for lam in partitions(M):
        N=sum(kostka(lam,m) for m in mus); rel=sum(1 for m in mus if kostka(lam,m)>0)
        tot+=N*flam(lam)
        print(lam, 'f=',flam(lam),'N=',N,'relorb=',rel,'f*rel=',flam(lam)*rel)
    print('orbits',len(O),'check sum f*N =',tot)
