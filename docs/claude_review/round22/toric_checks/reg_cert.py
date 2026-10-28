# Certify t_max(m) <= floor(10m/3) for k=5: Poly_{<=T}|_{S_m} has rank |S_m| mod p (=> over Q).
# t_max(m) = max ord_e f over Laurent polynomials f with Newton polytope in m*R_5 (sections of m(-K)).
# f vanishes to order >= t at e iff its coefficient vector is orthogonal to Poly_{<t}|_S; so
# full rank of Poly_{<=T}|_S  <=>  no f of order >= T+1.
# Memory-light: random subset of |S|+40 monomial columns of degree <= T; full rank of the subset suffices.
import itertools, sys, time, random
from flint import nmod_mat
p=(1<<61)-1
random.seed(1)
def S(m):
    return [x for x in itertools.product(range(-m,m+1),repeat=4) if sum(v for v in list(x)+[-sum(x)] if v>0)<=m]
for m in map(int,sys.argv[1:]):
    pts=S(m); N=len(pts); T=(10*m)//3
    ms=[a for a in itertools.product(range(T+1),repeat=4) if sum(a)<=T]
    cols=ms   # all monomials (a random subset of |S|+40 columns was insufficient at m=6: rank 5035)
    t0=time.time()
    M=nmod_mat(N,len(cols),p)
    for i,x in enumerate(pts):
        xs=[[pow(x[j]%p,e,p) for e in range(T+1)] for j in range(4)]
        for c,a in enumerate(cols):
            M[i,c]=xs[0][a[0]]*xs[1][a[1]]%p*xs[2][a[2]]%p*xs[3][a[3]]%p
    t1=time.time()
    r=M.rank()
    print(f"m={m} |S|={N} T={T} cols={len(cols)} of {len(ms)} rank={r} full={r==N} [build {t1-t0:.0f}s, rank {time.time()-t1:.0f}s]",flush=True)
