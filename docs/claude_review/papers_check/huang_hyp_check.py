# Check Huang's Hypothesis (**) (equivalent to (*): Eff simplicial, Lemma 6.3 p.35)
# for the permutohedral (Losev-Manin) fan in Z^k/Z1, k=3,4,5, over the splitting field.
import itertools
from fractions import Fraction
def solve(B, v):
    n=len(B); A=[[Fraction(B[j][i]) for j in range(n)]+[Fraction(v[i])] for i in range(n)]
    for c in range(n):
        p=next(r for r in range(c,n) if A[r][c]!=0); A[c],A[p]=A[p],A[c]
        for r in range(n):
            if r!=c and A[r][c]!=0:
                f=A[r][c]/A[c][c]; A[r]=[a-f*b for a,b in zip(A[r],A[c])]
    return [A[i][n]/A[i][i] for i in range(n)]
for k in (3,4,5):
    def vec(A): # coordinates in basis e_1..e_{k-1} of Z^k/Z1 (e_k = -sum)
        v=[0]*(k-1)
        for a in A:
            if a<k-1: v[a]+=1
            else: v=[x-1 for x in v]
        return v
    rays=[frozenset(A) for r in range(1,k) for A in itertools.combinations(range(k),r)]
    good=0; total=0
    for perm in itertools.permutations(range(k)):
        cone=[frozenset(perm[:i]) for i in range(1,k)]
        B=[vec(A) for A in cone]; total+=1
        ok=True
        for R in rays:
            if R in cone: continue
            c=solve(B,vec(R))
            if any(x>0 or x.denominator!=1 for x in c): ok=False; break
        good+=ok
    # centred primitive collections {A, A^c}: -K degree 2
    print(f"k={k}: rays={len(rays)}, geometric Pic rank={len(rays)-(k-1)}, maximal cones={total}, cones satisfying (**)={good}")
