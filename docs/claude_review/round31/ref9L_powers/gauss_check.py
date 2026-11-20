#!/usr/bin/env python3
"""Exact checks of the bounded-exponent Gaussian triangle (IX.0, powers paragraph).
Gaussian integers are pairs (x,y) of Python ints."""
import itertools, sys
from math import gcd, log
from sympy import factorint

def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def norm(a): return a[0]*a[0]+a[1]*a[1]
def gpow(a,k):
    r=(1,0)
    for _ in range(k): r=mul(r,a)
    return r
def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def divexact(a,b):
    # a/b in Z[i], assert exact
    n=norm(b); c=mul(a,conj(b))
    assert c[0]%n==0 and c[1]%n==0, ("not divisible",a,b)
    return (c[0]//n,c[1]//n)
def divides(b,a):
    n=norm(b); c=mul(a,conj(b)); return c[0]%n==0 and c[1]%n==0

def generator(p):
    for x in range(1,int(p**0.5)+1):
        y2=p-x*x; y=int(round(y2**0.5))
        if y*y==y2 and y>0: return (x,y)
    raise ValueError

def vpi(z,pi):
    k=0
    while divides(pi,z): z=divexact(z,pi); k+=1
    return k

def run(exps, eta=(1,0), verbose=True):
    primes=sorted(exps)
    pis={p:generator(p) for p in primes}
    for p in primes:  # pi and conj(pi) non-associate
        assert not divides(pis[p],conj(pis[p]))
    W=sum(e*log(p) for p,e in exps.items())
    allocs=list(itertools.product(*[range(exps[p]+1) for p in primes]))
    def point(a):
        z=eta
        for p,ap in zip(primes,a):
            z=mul(z,mul(gpow(pis[p],ap),gpow(conj(pis[p]),exps[p]-ap)))
        return z
    Npt=1
    for p,e in exps.items(): Npt*=p**e
    pts={a:point(a) for a in allocs}
    assert all(norm(z)==Npt for z in pts.values())
    assert len(set(pts.values()))==len(allocs)
    def H(a,b):
        h=(1,0)
        for k,p in enumerate(primes):
            h=mul(h,mul(gpow(pis[p],min(a[k],b[k])),gpow(conj(pis[p]),exps[p]-max(a[k],b[k]))))
        return h
    twoi_eta=mul((0,2),eta)
    # pair statement
    npairs=0; tvals={}
    for a,b in itertools.combinations(allocs,2):
        d=sub(pts[b],pts[a]); h=H(a,b)
        t=divexact(d,mul(twoi_eta,h))
        assert t[1]==0 and t[0]!=0, ("t not nonzero rational",a,b,t)
        # log|H| = (W-D)/2
        D=sum(abs(a[k]-b[k])*log(p) for k,p in enumerate(primes))
        assert abs(0.5*log(norm(h)) - (W-D)/2)<1e-9
        # the quotients are conjugate
        ua=divexact(pts[a],mul(eta,h)); ub=divexact(pts[b],mul(eta,h))
        assert ub==conj(ua)
        tvals[(a,b)]=t[0]; tvals[(b,a)]=-t[0]; npairs+=1
    # triangles
    ntri=0; nhidden_tri=0; ntwoblock=0
    for i,j,k in itertools.permutations(allocs,3):
        a1,a2,a3=i,j,k
        t12,t23,t13=tvals[(a1,a2)],tvals[(a2,a3)],tvals[(a1,a3)]
        H12,H23,H13=H(a1,a2),H(a2,a3),H(a1,a3)
        # raw identity H12 t12 + H23 t23 = H13 t13
        assert add(mul(H12,(t12,0)),mul(H23,(t23,0)))==mul(H13,(t13,0))
        # componentwise gcd G3
        G=(1,0)
        for q,p in enumerate(primes):
            m=min(a1[q],a2[q],a3[q]); M=max(a1[q],a2[q],a3[q])
            G=mul(G,mul(gpow(pis[p],m),gpow(conj(pis[p]),exps[p]-M)))
        A12,A23,A13=divexact(H12,G),divexact(H23,G),divexact(H13,G)
        assert add(mul(A12,(t12,0)),mul(A23,(t23,0)))==mul(A13,(t13,0))
        hidden=[]
        for q,p in enumerate(primes):
            x1,x2,x3=a1[q],a2[q],a3[q]
            pi=pis[p]; pb=conj(pi)
            ex={nm:(vpi(A,pi),vpi(A,pb)) for nm,A in (('12',A12),('23',A23),('13',A13))}
            if x1==x2!=x3:
                hidden.append(p)
                d=abs(x1-x3)
                exp12=(d,0) if x3<x1 else (0,d)
                assert ex['12']==exp12 and ex['23']==(0,0) and ex['13']==(0,0), (p,x1,x2,x3,ex)
            if len({x1,x2,x3})==3:
                nb=sum(1 for v in ex.values() if v!=(0,0));
                if nb==2: ntwoblock+=1
        if hidden:
            nhidden_tri+=1
            # hidden product X = part of A12 at hidden primes; recover from equation
            U=(1,0); X=(1,0)
            for q,p in enumerate(primes):
                pi=pis[p]; pb=conj(pi)
                fac=mul(gpow(pi,vpi(A12,pi)),gpow(pb,vpi(A12,pb)))
                if p in hidden: X=mul(X,fac)
                else: U=mul(U,fac)
            rhs=sub(mul(A13,(t13,0)),mul(A23,(t23,0)))
            Xrec=divexact(rhs,mul(U,(t12,0)))
            assert Xrec==X
            assert sorted(factorint(norm(Xrec)).keys())==sorted(hidden)
        ntri+=1
    if verbose:
        print(f"exps={exps} allocations={len(allocs)} pairs={npairs} ordered triples={ntri} "
              f"unordered={ntri//6} with hidden slot={nhidden_tri} (ordered) three-level slots in two blocks={ntwoblock}")
    return len(allocs),npairs,ntri//6

if __name__=="__main__":
    for eta in [(1,0),(0,1)]:
        run({5:2,13:2,17:2},eta)
    run({5:3,13:1,17:2})
    run({5:3,13:3})
    run({5:4,29:2})
    print("all assertions passed")
