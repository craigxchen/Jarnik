# Explicit integer 8-point cluster from the cube template at y=(-9/4,-2/3,1/4,8/3), t0 = p/2 in the optimal class;
# exact verification; best k-subclusters' asymptotic constants.
import sys, itertools, math
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round30/bridge_checks')
from localgcd import *
from cube_search import roots_points
from template import C_at
import mpmath as mp
mp.mp.dps=60
y=[Fr(-9,4),Fr(-2,3),Fr(1,4),Fr(8,3)]
pr=roots_points(y,[1,-1,-1,1]); zs=[prod_lin(rr) for rr in pr]
target=None
print("C_inf (8 points):",C_inf(pr,8)[0])
# find p in the optimal class with |p| ~ 10^4 ... 10^30
for p0 in (10**4, 10**12, 10**30):
    best=None
    for p in range(p0, p0+4000):
        if p%2==0: continue
        C,gabs,g,den=C_at(zs,Fr(p,2))
        if best is None or C<best[0]: best=(C,p,g,den)
    C,p,g,den=best
    t0=Fr(p,2); vals=[peval(z,t0) for z in zs]
    # normalised Gaussian integers w = val/(g/den)
    gg=G(Fr(g[0],den),Fr(g[1],den))
    W=[v/gg for v in vals]
    assert all(w.a.denominator==1 and w.b.denominator==1 for w in W)
    Wi=[(int(w.a),int(w.b)) for w in W]
    n=Wi[0][0]**2+Wi[0][1]**2
    assert all(a*a+b*b==n for a,b in Wi) and len(set(Wi))==8
    h=(0,0)
    for w in Wi: h=gi_gcd(h,w) if h!=(0,0) else w
    assert h[0]**2+h[1]**2==1
    print(f"t0={p}/2: R^2={n} (digits {len(str(n))}), C(t0)={mp.nstr(C,15)}")
    if p0==10**4:
        print("  points:",Wi); print("  factor R^2:",sp.factorint(n))
# sub-clusters
for k in (5,6,7):
    bestk=None
    for S in itertools.combinations(range(8),k):
        out=C_inf([pr[i] for i in S],8)
        if out and (bestk is None or out[0]<bestk[0]): bestk=(out[0],S)
    print(f"best {k}-subcluster C_inf = {bestk[0]:.6g} subset {bestk[1]}",flush=True)
