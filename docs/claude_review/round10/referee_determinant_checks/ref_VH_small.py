# Referee: at the author's small evidence sizes, is V_H (span of integer kernel vectors with |c|_inf <= H,
# H = (U T^{m-1})^2) equal to the whole V(A,m)?  If an integer basis of V(A,m) has sup-norm <= H, yes.
# Same slab-box construction as determinant_checks/baselocus2.py (seed 1), exact rational nullspace.
import numpy as np, itertools, math, sys
from fractions import Fraction
import sympy as sp
rng = np.random.default_rng(1)
def slab_box(n,T,r,frac=0.5):
    d=n-1
    box=np.array(list(itertools.product(range(-T,T+1),repeat=d)))
    beta=rng.choice([-1,1],size=(r,n)).astype(float)
    X=box@(beta[:,:d]-beta[:,[d]]).T
    mx=np.abs(X).max(axis=1); w=np.sort(mx)[int(frac*len(mx))]
    return box[mx<=w].astype(int), beta, w
out=[]
for (n,T,r,m) in [(4,3,8,5),(4,4,8,6),(5,2,10,5)]:
    A,beta,w=slab_box(n,T,r)
    d=n-1; U=len(A)
    mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=d) if sum(a)==t]
    M=sp.Matrix([[int(np.prod([int(k[j])**a[j] for j in range(d)])) for k in A] for a in mons])
    ns=M.nullspace()
    mx=0
    for v in ns:
        den=sp.ilcm(*[sp.fraction(x)[1] for x in v]); c=[int(x*den) for x in v]
        g=0
        for x in c: g=math.gcd(g,abs(x))
        c=[x//g for x in c]
        mx=max(mx,max(abs(x) for x in c))
    H=(U*T**(m-1))**2
    s=f"n={n} T={T} m={m} |A|={U} dimV={len(ns)} max sup-norm of echelon integer basis={mx} H={H} V_H=V(A,m)? {mx<=H}"
    print(s,flush=True); out.append(s)
open(sys.argv[1],"w").write("\n".join(out)+"\n")
