# Template constant of a Q(i)(t)-cluster: C(t0) = maxchord/sqrt(|z(t0)| |g(t0)|), g = fractional Gaussian gcd of values.
# Asymptotically C -> max|gamma| / sqrt|g(t0)|.  Scan t0 = p/q.
import sys, itertools, math
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
from glib import *
from cube_basic import cube
import mpmath as mp
mp.mp.dps=50

def values(pts,t0):
    return [peval(z,t0) for z in pts]

def C_at(pts,t0):
    v=values(pts,t0)
    g,den=frac_gcd(v)
    gabs=mp.sqrt(mp.mpf(gnorm_int(g)))/den
    R=mp.sqrt(mp.mpf(v[0].norm().numerator)/v[0].norm().denominator)
    ch=0
    for a,b in itertools.combinations(v,2):
        d=a-b; n=d.norm(); c=mp.sqrt(mp.mpf(n.numerator)/n.denominator)
        ch=max(ch,c)
    # arc length rather than chord: arc = 2R asin(ch/(2R))
    arc=2*R*mp.asin(ch/(2*R))
    return arc/mp.sqrt(R*gabs), gabs, g, den

def gammas(pts,D):
    out=[]
    for a,b in itertools.combinations(pts,2):
        d=psub(a,b); out.append(d[D//2] if len(d)>D//2 else G(0))
    return out

if __name__=="__main__":
    y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
    roots,P=cube(y); pts=[P[k][0] for k in P]
    gm=max(mp.sqrt(mp.mpf(g.norm().numerator)/g.norm().denominator) for g in gammas(pts,8))
    print("max|gamma| =",gm)
    best=None
    for q in [1,2,3,4,6,8,9,12,18,24,36,72]:
        bq=None
        for p in range(10**7, 10**7+2*72*q):
            if math.gcd(p,q)!=1: continue
            C,gabs,g,den=C_at(pts,Fr(p,q))
            if bq is None or C<bq[0]: bq=(C,p,q,gabs)
        print("q=",q," best C=",mp.nstr(bq[0],12)," at p mod ",bq[1]," |g|=",mp.nstr(bq[3],8), " gm/sqrt|g|=",mp.nstr(gm/mp.sqrt(bq[3]),12))
        if best is None or bq[0]<best[0]: best=bq
    print("BEST",best)
