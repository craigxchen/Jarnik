"""Independent CRT-lattice upper bound (referee).  k=8: labels 1..7 (anchor 0 at infinity), X_1 = 0.
Blocks: distinct random primes in [A,2A] on the 119 visible T.  Lattice of X=(X_2..X_7) with
X_i = X_j mod n_T for i,j in T.  Basis built by successive sublattice intersection, LLL by flint,
search coefficients in [-2,2]^6.  Report min H = sum_{i<j} log|t_ij|, verify divisibility exactly."""
import sys, random
from itertools import combinations, product
from math import log
import flint
from sympy import primerange
k=int(sys.argv[1]); A=int(sys.argv[2]); seed=int(sys.argv[3]); m=k-1
rng=random.Random(seed)
Ts=[frozenset(S) for s in range(2,m) for S in combinations(range(1,m+1),s)]
ps=rng.sample(list(primerange(A,2*A)),len(Ts)); n=dict(zip(Ts,ps))
dim=m-1   # coordinates X_2..X_m ; X_1=0
B=[[1 if i==j else 0 for j in range(dim)] for i in range(dim)]
def coord(vec,i): return 0 if i==1 else vec[i-2]
for T in Ts:
    q=n[T]; L=sorted(T)
    for a in L[1:]:
        # constraint X_a - X_{L0} = 0 mod q
        vals=[(coord(b,a)-coord(b,L[0]))%q for b in B]
        piv=next((i for i,v in enumerate(vals) if v),None)
        if piv is None: continue
        inv=pow(vals[piv],-1,q)
        for i in range(dim):
            if i!=piv and vals[i]:
                c=vals[i]*inv%q
                B[i]=[x-c*y for x,y in zip(B[i],B[piv])]
        B[piv]=[q*y for y in B[piv]]
M=flint.fmpz_mat(B)
idx=abs(int(M.det())); expect=1
for T in Ts: expect*=n[T]**(len(T)-1)
assert idx==expect, "index mismatch"
R=M.lll()
basis=[[int(R[i,j]) for j in range(dim)] for i in range(dim)]
best=None
for c in product(range(-2,3),repeat=dim):
    if not any(c): continue
    Y=[sum(c[l]*basis[l][r] for l in range(dim)) for r in range(dim)]
    X=[0]+Y
    if len(set(X))<m: continue
    H=0.0
    for i,j in combinations(range(1,m+1),2):
        G=1
        for T in Ts:
            if i in T and j in T: G*=n[T]
        d=X[i-1]-X[j-1]
        assert d%G==0
        H+=log(abs(d//G))
    if best is None or H<best[0]: best=(H,X)
w=log(A); eta=max(abs(log(p)-w) for p in ps)/w
print("k=%d A=%d seed=%d: index ok; min H=%.2f w=%.3f H/w=%.1f eta=%.3f max log10|X|=%.1f"%(k,A,seed,best[0],w,best[0]/w,eta,max(log(abs(x),10) for x in best[1] if x)),flush=True)
