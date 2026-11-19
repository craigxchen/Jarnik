# Independent reproduction of the finite Gini certificate, using the CDF formula
# GMD(mu) = 2 * sum over unit gaps of F(t)(1-F(t)) (integer-supported laws), exact rationals.
from fractions import Fraction as Fr
from itertools import product
def law(c, nonconst):
    d={}
    pats=[b for b in product((0,1),repeat=5) if (not nonconst) or (0<sum(b)<5)]
    for b in pats:
        x=sum(bi*ci for bi,ci in zip(b,c)); d[x]=d.get(x,0)+Fr(1,len(pats))
    return d
def gmd(d):
    lo,hi=min(d),max(d); F=Fr(0); s=Fr(0)
    for t in range(lo,hi):
        F+=d.get(t,0); s+=F*(1-F)
    return 2*s
def mix(d1,d2,r):
    keys=set(d1)|set(d2); return {k:(1-r)*d1.get(k,0)+r*d2.get(k,0) for k in keys}
def maxr(c):
    U=law(c,False); Us=law(c,True)
    # GMD(mix) is a quadratic in r (not obviously; check by fitting 3 points and verifying 2 more)
    g=lambda r: gmd(mix(U,Us,r))
    g0,gh,g1=g(Fr(0)),g(Fr(1,2)),g(Fr(1))
    # quadratic through these: A r^2 + B r + C
    C=g0; A=2*g1+2*g0-4*gh; B=g1-g0-A
    for r in (Fr(1,3),Fr(3,7),Fr(5,6)):
        assert g(r)==A*r*r+B*r+C, "not quadratic"
    vals=[g0,g1]
    if A<0:
        rs=-B/(2*A)
        if 0<rs<1: vals.append(g(rs))
    return max(vals)
def cent(x,q):
    x%=q
    return x-q if x>q//2 else x
for q in (3,5,7):
    best=Fr(0); worstc=None; ncls=set()
    for c in product(range(q),repeat=5):
        if not any(c) or sum(c)%q: continue
        key=min(tuple(sorted(cent(k*x,q) for x in c)) for k in range(1,q))
        if key in ncls: continue
        ncls.add(key)
        val=min(maxr([cent(k*x,q) for x in c]) for k in range(1,q))
        if val>best: best=val; worstc=c
    print(q,"classes",len(ncls),"worst min-over-scaling max-over-r GMD =",best,float(best),"q/2=",q/2,"worst c",worstc)
