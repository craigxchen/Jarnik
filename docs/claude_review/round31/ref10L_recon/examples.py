from math import gcd, isqrt
from itertools import combinations, permutations

def H(a,b,c): return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)

def check_table(primes, sig, d, delta, C2):
    """sig[i]: tuple of bits; d[(i,j)] magnitude; delta[p]; C2 = C^2 (rational as int pair ok)"""
    q=1
    for p in primes: q*=p
    ok=True
    n=len(sig)
    for i,j in combinations(range(n),2):
        x=d[(i,j)]
        h=1
        for k,p in enumerate(primes):
            if sig[i][k]==sig[j][k]: h*=p
            else:
                if x%p not in (delta[p]%p,(-delta[p])%p): ok=False; print("res fail",i,j,p)
        if gcd(x,q)!=h: ok=False; print("gcd fail",i,j,x,h)
        if x*x> C2*q: ok=False; print("cap fail",i,j)
    return ok,q

def line_embeddable(n,d):
    # brute force: fix vertex 0 at 0, try signs
    for signs in range(1<<(n-1)):
        y=[0]+[ (d[(0,i)] if (signs>>(i-1))&1 else -d[(0,i)]) for i in range(1,n)]
        if all(abs(y[i]-y[j])==d[(i,j)] for i,j in combinations(range(n),2)): return y
    return None

# Example 1: q=17*23*29, 4-cycle
primes=[17,23,29]
D={(0,1):17,(1,2):29,(2,3):17,(0,3):29,(0,2):46,(1,3):46}
sig=[(0,0,0),(0,1,1),(1,0,1),(1,1,0)]
# at 17: {0,1},{2,3}; at 23:{0,2},{1,3}; at 29: {1,2},{0,3}
delta={17:12,23:6,29:17}
ok,q=check_table(primes,sig,D,delta,1)
print("ex1 ok",ok,"q",q,"isqrt",isqrt(q))
for t in combinations(range(4),3):
    a,b,c=D[(t[0],t[1])],D[(t[0],t[2])],D[(t[1],t[2])]
    print(t,(a,b,c),H(a,b,c))
print("line embeddable:",line_embeddable(4,D))

# Example 2
a,b=83,59
u=a*(a+b); v=b*(a+2*b)
Y=[0,u,v,u+v,(a+b)*(a+2*b)]
primes=[59,67,71,83]; q=59*67*71*83
print("q",q, "marks",Y)
for p in primes: print(p, sorted(set(y%p for y in Y)), [y%p for y in Y])
sig=[]; delta={}
for p in primes:
    rs=sorted(set(y%p for y in Y)); assert len(rs)==2
    delta[p]=(rs[1]-rs[0])%p
sig=[tuple(0 if y%p==sorted(set(z%p for z in Y))[0] else 1 for p in primes) for y in Y]
print("sigs",sig, "distinct", len(set(sig))==5)
D={(i,j):abs(Y[i]-Y[j]) for i,j in combinations(range(5),2)}
print("genuine table ok (C=6):",check_table(primes,sig,D,delta,36)[0])
for e,x in D.items(): print(e,x,gcd(x,q))
D2=dict(D); D2[(0,3)]=abs(u-v)
print("modified edge",D2[(0,3)],"ok:",check_table(primes,sig,D2,delta,36)[0], "max^2/q", max(D2.values())**2/q)
for t in combinations(range(5),3):
    a_,b_,c_=D2[(t[0],t[1])],D2[(t[0],t[2])],D2[(t[1],t[2])]
    hh=H(a_,b_,c_)
    print(t,(a_,b_,c_),hh, hh%(q*q)==0 and hh//(q*q))
print("line embeddable modified:",line_embeddable(5,D2))
# max embeddable subsets
for r in (4,):
    for S in combinations(range(5),r):
        idx={k:i for i,k in enumerate(S)}
        dd={(idx[i],idx[j]):D2[(i,j)] for i,j in combinations(S,2)}
        print(S, line_embeddable(r,dd) is not None)
