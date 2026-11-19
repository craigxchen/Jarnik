# brute force: all points A+B eta with max(deg A, deg B+1) = D (all leading coefficients), group by exact norm,
# then max clique in the graph "pole order of difference <= floor(D/2)" by exhaustive search (no ultrametric shortcut).
import itertools, sys
from collections import defaultdict
P=3
def mul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%P
    return r
def tr(a):
    a=[x%P for x in a]
    while a and a[-1]==0: a.pop()
    return tuple(a)
def deg(a): return len(tr(a))-1
Delta=[1,0,2]
def sqf(N):
    import numpy as np
    from fp import trim, squarefree
    return squarefree(trim(np.array(N)))
def maxclique(pts, ok):
    best=1 if pts else 0
    n=len(pts)
    def rec(cl, cand):
        nonlocal best
        if len(cl)>best: best=len(cl)
        for i,c in enumerate(cand):
            if len(cl)+len(cand)-i<=best: return
            rec(cl+[c],[d for d in cand[i+1:] if ok(c,d)])
    rec([], list(range(n)))
    return best
for D in range(1,int(sys.argv[1])+1):
    H=D//2
    fib=defaultdict(list)
    for A in itertools.product(range(3),repeat=D+1):
        for B in itertools.product(range(3),repeat=D):
            dA, dB = deg(A), deg(B)
            if max(dA, (dB+1) if dB>=0 else -1)!=D: continue
            N=tr([x-y for x,y in itertools.zip_longest(mul(A,A),mul(Delta,mul(B,B)),fillvalue=0)])
            fib[N].append((A,B))
    def ok(i,j,pts):
        (A1,B1),(A2,B2)=pts[i],pts[j]
        dA=deg([x-y for x,y in zip(A1,A2)]); dB=deg([x-y for x,y in zip(B1,B2)])
        po=max(dA, dB+1 if dB>=0 else -1)
        return po<=H
    ba=bs=0
    for N,pts in fib.items():
        assert len(N)-1==2*D
        if len(pts)<=ba and len(pts)<=bs: continue
        c=maxclique(pts, lambda i,j: ok(i,j,pts))
        ba=max(ba,c)
        if c>bs and sqf(list(N)): bs=c
    print(D, 'all', ba, 'sf', bs, flush=True)
