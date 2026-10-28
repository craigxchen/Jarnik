# Independent check of isotypic maximal orders t_+(m), t_-(m) on S_m = m R_5 ∩ M.
# ord_e(sum c_s x^s) >= t  iff  sum_s c_s p(s) = 0 for all polys p of degree < t (Lemma 3.1 of toric.md).
# iota-even f (c_s = c_{-s}): only even-degree p matter; iota-odd: only odd-degree p.
# Upper bound certified by full rank mod prime P (full rank mod P => full rank over Q).
import sys, itertools
from flint import nmod_mat
P = (1<<61)-1
def S(m):
    pts=[]
    for s in itertools.product(range(-m,m+1),repeat=4):
        s5=-sum(s); v=s+(s5,)
        if sum(abs(a) for a in v)<=2*m: pts.append(v)
    return pts
def monos(d,parity):
    out=[]
    for k in range(d+1):
        if k%2!=parity: continue
        for e in itertools.combinations_with_replacement(range(4),k):
            out.append(e)
    return out
def rank(reps,d,parity):
    rows=monos(d,parity)
    M=nmod_mat(len(rows),len(reps),P)
    for i,e in enumerate(rows):
        for j,s in enumerate(reps):
            v=1
            for a in e: v=v*s[a]
            M[i,j]=v%P
    return M.rank()
def tmax(reps,parity,guess):
    # smallest d with full rank over parity-d monomials
    n=len(reps)
    d=guess-1
    r_lo=rank(reps,d,parity) if d>=0 else 0
    r_hi=rank(reps,guess,parity)
    return r_lo,r_hi,n
for m in map(int,sys.argv[1:]):
    pts=S(m); st=set(pts)
    reps=[]; seen=set()
    for s in pts:
        if s in seen: continue
        seen.add(s); seen.add(tuple(-a for a in s)); reps.append(s)
    zero=(0,)*5
    even=reps                      # basis x^s+x^{-s} (and 1)
    odd=[s for s in reps if s!=zero]  # basis x^s-x^{-s}
    t=(10*m)//3
    out=[f"m={m} |S|={len(pts)} n_even={len(even)} n_odd={len(odd)} floor(10m/3)={t}"]
    for name,reps_,par in (("even",even,0),("odd",odd,1)):
        # claimed: part with parity of t reaches t, the other t-1
        claim = t if t%2==par else t-1
        rlo,rhi,n=tmax(reps_,par,claim)
        out.append(f"  {name}: claim t={claim}: rank(deg<={claim-1})={rlo} (deficient={rlo<n}), rank(deg<={claim})={rhi} (full={rhi==n}) n={n}")
    print("\n".join(out),flush=True)
