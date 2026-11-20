"""Tests at moduli where genuine 5-mark configurations exist (projective family a,b).
For q=product of primes dividing the marks' structure, with delta read off from the genuine
configuration, compute B_free (max clique of Cayley graph of admissible patterns),
B_int (best eps), all admissible magnitudes, and largest line-embeddable subsets of
random tables on maximum cliques."""
import random, sys
from math import gcd, isqrt
from itertools import product, combinations
sys.path.insert(0,'.')
from stmt1 import max_clique
from tables import embed, max_embeddable

def crt(rs, ms):
    x,M=0,1
    for r,m in zip(rs,ms):
        t=((r-x)*pow(M,-1,m))%m
        x+=M*t; M*=m
    return x%M

def admissible(primes, deltas, C2):
    k=len(primes); q=1
    for p in primes:q*=p
    L=isqrt(C2*q)
    al={}
    for pat in range(1,1<<k):
        S=[t for t in range(k) if (pat>>t)&1]
        A=[t for t in range(k) if not (pat>>t)&1]
        lst=set()
        for sg in product((1,-1),repeat=len(S)):
            rs=[0]*len(A)+[(sg[i]*deltas[t])%primes[t] for i,t in enumerate(S)]
            ms=[primes[t] for t in A]+[primes[t] for t in S]
            x=crt(rs,ms)
            for y in (x,):
                if 1<=y<=L: lst.add(y)
        al[pat]=sorted(lst)
    return al,q,L

def Bint(primes,deltas,L):
    best=0
    for eps in product((1,-1),repeat=len(primes)):
        c=sum(1 for y in range(L+1) if all(y%p in (0,(e*d)%p) for p,e,d in zip(primes,eps,deltas)))
        best=max(best,c)
    return best

def study(primes,deltas,C2,trials=20,label=""):
    al,q,L=admissible(primes,deltas,C2)
    k=len(primes); n=1<<k
    adj=[[i!=j and bool(al[i^j]) for j in range(n)] for i in range(n)]
    cl=max_clique(n,adj)
    bi=Bint(primes,deltas,L)
    multi={pat:l for pat,l in al.items() if len(l)>1}
    hs={}
    for pat in multi:
        h=1
        for t,p in enumerate(primes):
            if not (pat>>t)&1: h*=p
        hs[pat]=h
    res={}
    for tr in range(trials):
        ch={}
        def d(i,j):
            key=(min(i,j),max(i,j))
            if key not in ch: ch[key]=random.choice(al[i^j])
            return ch[key]
        r,T=max_embeddable(cl,d)
        res[len(cl)-r]=res.get(len(cl)-r,0)+1
    print(label,"q",q,"C2",C2,"L",L,"B_free",len(cl),"B_int",bi,
          "ok" if len(cl)<=2*max(4,bi) else "VIOLATION",
          "patterns with >1 magnitude: h values",sorted(set(hs.values())),
          "deletions needed (random tables):",res,flush=True)
    return len(cl),bi,cl,al

random.seed(3)
a,b=83,59
primes=[59,67,71,83]
Y=[0,a*(a+b),b*(a+2*b),a*(a+b)+b*(a+2*b),(a+b)*(a+2*b)]
deltas=[]
for p in primes:
    rs=sorted(set(y%p for y in Y)); deltas.append(rs[1]-rs[0])
for C2 in (36,49,58):
    study(primes,deltas,C2,label="5-mark modulus")
# random deltas, same primes
for C2 in (36,58):
    mx=(0,0)
    for _ in range(40):
        dl=[random.randint(1,(p-1)//2) for p in primes]
        bf,bi,_,_=study(primes,dl,C2,trials=3,label="random delta")
print("clique at C2=36 for 5-mark modulus:")
bf,bi,cl,al=study(primes,deltas,36,trials=1,label="again")
print([format(x,'04b') for x in cl])
for i,j in combinations(cl,2): print(format(i,'04b'),format(j,'04b'),al[i^j])
