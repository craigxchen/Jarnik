import itertools, sympy as sp
def fake(k,c):
    D=2**k; L=D//2-c-1
    us=list(itertools.product([0,1],repeat=k))
    nz=[u for u in us if any(u)]
    # choose B: basis vectors first then others, |B|=L
    basis=[tuple(1 if i==j else 0 for i in range(k)) for j in range(k)]
    rest=[u for u in nz if u not in basis]
    B=basis+rest[:L-k]
    assert len(B)==L
    chi=lambda b,u: (-1)**sum(x*y for x,y in zip(b,u))
    W=sp.Matrix([[chi(b,u) for u in us] for b in B])
    cols=[tuple(W[:,j]) for j in range(D)]
    distinct=len(set(cols))==D
    pts=[s for s in itertools.product([0,1],repeat=D) if all(x==0 for x in W*sp.Matrix(s))]
    md=min(sum(a!=b for a,b in zip(s,t)) for s,t in itertools.combinations(pts,2))
    return D,L,len(pts),D+2*c+2,md,distinct
for k,c in [(3,0),(4,0),(4,1),(4,2)]:
    print(fake(k,c))
