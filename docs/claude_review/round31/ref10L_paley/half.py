# Half-row set R = {i : K_i0 = 1} (= nonresidues). Find zero-sum a on R orthogonal to many columns S_v|R.
import sys, random, math
from itertools import combinations
from fractions import Fraction
from paley1 import chi, S_of
from spark import rank_mod
def nullvec(rows,n):
    # exact rational null vector of a (n-1) x n system of rank n-1, return primitive integer vector
    A=[[Fraction(x) for x in r] for r in rows]; piv=[];r=0
    for c in range(n):
        p=next((i for i in range(r,len(A)) if A[i][c]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]
        A[r]=[x/A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][c]!=0:
                f=A[i][c]; A[i]=[A[i][t]-f*A[r][t] for t in range(n)]
        piv.append(c); r+=1
    if r!=n-1: return None
    free=[c for c in range(n) if c not in piv][0]
    x=[Fraction(0)]*n; x[free]=Fraction(1)
    for i,c in enumerate(piv): x[c]=-A[i][free]
    L=1
    for t in x: L=L*t.denominator//math.gcd(L,t.denominator)
    xi=[int(t*L) for t in x]; g=0
    for t in xi: g=math.gcd(g,t)
    return tuple(t//g for t in xi)
def explore(m,exh,samples=20000,seed=1):
    random.seed(seed)
    S=S_of(m); R=[i for i in range(m) if chi(0-i,m)==1]; n=len(R)
    cols=[tuple(S[i][v] for i in R) for v in range(m)]
    one=tuple([1]*n)
    found={}
    it=combinations(range(m),n-2) if exh else (random.sample(range(m),n-2) for _ in range(samples))
    for Z0 in it:
        rows=[one]+[cols[v] for v in Z0]
        a=nullvec(rows,n)
        if a is None: continue
        if a[next(i for i in range(n) if a[i])]<0: a=tuple(-t for t in a)
        if a in found: continue
        Z=[v for v in range(m) if sum(a[i]*cols[v][i] for i in range(n))==0]
        found[a]=len(Z)
    best=max(found.values())
    hist={}
    for z in found.values(): hist[z]=hist.get(z,0)+1
    # do annihilators with |Z| > m/2 span the zero-sum space (dim n-1)?
    good=[a for a,z in found.items() if z>m/2]
    rk=rank_mod([list(a) for a in good]) if good else 0
    l1=min((sum(abs(t) for t in a) for a in good),default=None)
    print("m=%d |R|=%d: max|Z|=%d (m/2=%.1f, generic (m-5)/2=%d); hist %s; #a with |Z|>m/2: %d, their rank %d of %d; min l1 %s"%(m,n,best,m/2,(m-5)//2,dict(sorted(hist.items())),len(good),rk,n-1,l1))
explore(11,True); explore(19,True); explore(23,True)
explore(31,False,20000); explore(43,False,5000)
