# Quadratic profile: bits x_i, x_i x_j (i<j), and all nonconstant linear parities u.x. Check Sidon and real affine rank.
from itertools import combinations
from sympy import Matrix
def dot(u,x): return bin(u&x).count('1')&1
for r in (3,4,5,6):
    n=1<<r
    cols=[lambda x,i=i:(x>>i)&1 for i in range(r)]
    cols+=[lambda x,i=i,j=j:((x>>i)&(x>>j))&1 for i,j in combinations(range(r),2)]
    cols+=[lambda x,u=u:dot(u,x) for u in range(1,n)]
    s=len(cols)
    sig=[sum(c(x)<<k for k,c in enumerate(cols)) for x in range(n)]
    sums={}
    ok=True
    for x,y in combinations(range(n),2):
        h=sig[x]^sig[y]
        if h in sums: ok=False
        sums[h]=(x,y)
    D=len(sums)+1
    rows=[[1]+[c(x) for c in cols] for x in range(n)]
    rk=Matrix(rows).rank() if r<=5 else 'skip'
    # F2 rank of signature matrix (with constant)
    vs=[(sig[x]<<1)|1 for x in range(n)]; basis=[]
    for v in vs:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    print(f"r={r} n={n} s={s} Sidon={ok} D={D} (1+C(n,2)={1+n*(n-1)//2}) real rank(with const)={rk} F2 rank={len(basis)} (bound 1+r+C(r,2)={1+r+r*(r-1)//2})")
