"""Referee's independent library (no import of author code).
Points of P^1(Q) as primitive integer vectors; labels 0..k-1.
Splits: frozenset S of labels with 0 not in S, 2<=|S|<=k-2."""
from itertools import combinations
from math import gcd, log
from sympy import factorint

def splits(k):
    out=[]
    for s in range(2,k-1):
        for S in combinations(range(1,k),s):
            out.append(frozenset(S))
    return out

def det(u,v): return u[0]*v[1]-u[1]*v[0]

def vp(n,p):
    n=abs(n); c=0
    while n%p==0: n//=p; c+=1
    return c

def contacts_minformula(P):
    """Thm B(iv) min-formula, by factoring all minors (small inputs only)."""
    k=len(P); primes=set()
    for i,j in combinations(range(k),2):
        d=det(P[i],P[j]); assert d!=0
        primes|=set(factorint(abs(d)).keys())
    Sp=splits(k); out={S:1 for S in Sp}; lp={}
    for p in primes:
        g={}
        for i,j in combinations(range(k),2):
            g[(i,j)]=g[(j,i)]=vp(det(P[i],P[j]),p)
        for S in Sp:
            C=[c for c in range(k) if c not in S]
            best=None
            for a,b in combinations(sorted(S),2):
                for c,d in combinations(C,2):
                    v=min(g[a,b]+g[c,d]-g[a,c]-g[b,d], g[a,b]+g[c,d]-g[a,d]-g[b,c])
                    best=v if best is None else min(best,v)
            l=max(0,best)
            if l: out[S]*=p**l; lp[(S,p)]=l
    return out, lp

def contacts_classes(x):
    """integers x_0..x_{k-1} (points (x_i,1); infinity unmarked).  Residue-class rule:
    l_split(p) = #{t>=1: S is a class mod p^t} + #{t>=1: S^c is a class mod p^t}."""
    k=len(x); Sp=set(splits(k)); out={S:1 for S in Sp}
    primes=set()
    for i,j in combinations(range(k),2):
        primes|=set(factorint(abs(x[i]-x[j])).keys())
    for p in primes:
        t=1
        while True:
            q=p**t; cl={}
            for i,xi in enumerate(x): cl.setdefault(xi%q,[]).append(i)
            big=[frozenset(c) for c in cl.values() if len(c)>=2]
            if not big: break
            for C in big:
                if len(C)<=k-2:
                    S=C if 0 not in C else frozenset(range(k))-C
                    out[S]*=p
            t+=1
    return out

def contacts_gcd_anchor(X):
    """anchor at infinity (label 0), X_1..X_m integers (labels 1..m).  g_T = gcd_{i,j in T}(X_i-X_j)
    = delta * prod_{T' >= T visible} d_T'; invert top-down. returns d (dict over visible T subset [m]), delta."""
    m=len(X); delta=0
    for i,j in combinations(range(m),2): delta=gcd(delta,X[i]-X[j])
    Ts=[frozenset(S) for s in range(m-1,1,-1) for S in combinations(range(1,m+1),s)]
    d={}
    for T in Ts:
        g=0; L=sorted(T)
        for a,b in combinations(L,2): g=gcd(g,X[a-1]-X[b-1])
        prod=delta
        for T2,v in d.items():
            if T<T2: prod*=v
        assert g%prod==0
        d[T]=g//prod
    return d,delta
