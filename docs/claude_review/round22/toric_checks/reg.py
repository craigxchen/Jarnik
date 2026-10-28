# reg(S_m) for S_m = mR cap M (R = root polytope of A_{k-1}), via ranks mod a large prime.
# reg(S) = min d with Poly_{<=d}|_S = F^S;  t_max(m) = max ord_e f over f in H^0(m(-K)) equals reg(S_m).
# rank mod p <= rank over Q, so "full rank mod p" certifies (over Q) that no f of order > d exists.
import itertools, sys, time
from flint import nmod_mat
p = (1<<61)-1  # prime
def S(k,m):
    pts=[]
    for x in itertools.product(range(-m,m+1),repeat=k-1):
        y=list(x)+[-sum(x)]
        if sum(v for v in y if v>0)<=m: pts.append(x)
    return pts
def monos(n,d):
    return [a for t in range(d+1) for a in itertools.product(range(t+1),repeat=n) if sum(a)==t] if False else \
           [a for a in itertools.product(range(d+1),repeat=n) if sum(a)<=d]
def rank_upto(pts,d):
    n=len(pts[0]); ms=monos(n,d)
    rows=[]
    for x in pts:
        rows.append([ (pow(x[0]%p,a[0],p)*pow(x[1]%p,a[1],p)%p*pow(x[2]%p,a[2],p)%p*(pow(x[3]%p,a[3],p) if n>3 else 1))%p for a in ms])
    M=nmod_mat(rows,p)
    return M.rank(), len(ms)
k=int(sys.argv[1]); mmax=int(sys.argv[2]); dstart=sys.argv[3] if len(sys.argv)>3 else None
for m in range(1,mmax+1):
    pts=S(k,m); N=len(pts)
    d=int(3*m) if k==5 else int(3*m)-1
    d=max(d-1,0)
    t0=time.time()
    while True:
        r,c=rank_upto(pts,d)
        if r==N: break
        d+=1
    # d = reg: rank full at d, deficient at d-1
    print(f"k={k} m={m} |S|={N} reg={d} reg/m={d/m:.4f}  (Vandermonde-type lower {('10/3' if k==5 else '3')}*m={ (10*m/3 if k==5 else 3*m):.3f}, mean-width upper {(15*m/4 if k==5 else 7*m/2):.3f})  [{time.time()-t0:.1f}s]",flush=True)
