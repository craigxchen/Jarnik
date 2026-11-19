"""Actual circles: 5-point windows (consecutive by angle), gcd-normalised.
Checks (4) on squareclass-dependent windows, checks (1) on all pairs, records index distribution."""
import sys, random, math
from itertools import combinations
from collections import Counter
from math import isqrt, atan2
from mpmath import mp, mpf, sqrt as msqrt
from sympy import primerange
sys.path.insert(0, sys.argv[0].rsplit('/', 1)[0])
from gauss_util import *

mp.dps = 50
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
P1 = [p for p in primerange(5, 400) if p % 4 == 1]


def points_on(N):
    pts = []
    for x in range(-isqrt(N), isqrt(N)+1):
        y2 = N - x*x
        y = isqrt(y2)
        if y*y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    pts = list(set(pts))
    pts.sort(key=lambda t: atan2(t[1], t[0]))
    return pts


idx = Counter()
par = Counter()
minratio_dep = None
minratio_all = None
nwin = 0
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 150):
    k = random.randint(3, 6)
    ps = random.sample(P1, k)
    N = 1
    for p in ps:
        N *= p**random.randint(1, 3)
    if N > 4*10**11:
        continue
    pts = points_on(N)
    n = len(pts)
    for s in range(n):
        z = [pts[(s+t) % n] for t in range(5)]
        g = z[0]
        for w in z[1:]:
            g = ggcd(g, w)
        zz = [gdivmod_exact(w, g) for w in z]
        H = gnorm(zz[0])
        if H == 1:
            continue
        prof, units = profile(zz, H)
        # check primitivity: min 0, max e
        for p, (e, a) in prof.items():
            assert min(a) == 0 and max(a) == e
        R = msqrt(mpf(H))
        Delta = max(msqrt(mpf((u[0]-v[0])**2+(u[1]-v[1])**2)) for u, v in combinations(zz, 2))/R
        # Liouville (1) on every pair
        for i, j in combinations(range(5), 2):
            dij = sum(abs(a[i]-a[j])*math.log(p) for p, (e, a) in prof.items())
            dist = msqrt(mpf((zz[i][0]-zz[j][0])**2+(zz[i][1]-zz[j][1])**2))/R
            assert dist >= mp.e**(-mpf(dij)/2)*(1-mpf(10)**-30), (zz, i, j)
        cols = columns(prof)
        I = lattice_index(cols)
        idx[I] += 1
        dep = parity_dependent(prof, units)
        dep1 = parity_dependent(prof, units, with_units=False)
        oneclass = len(set(units)) == 1
        par[(dep, dep1, oneclass)] += 1
        ratio = Delta*mpf(H)**0.25
        nwin += 1
        if minratio_all is None or ratio < minratio_all:
            minratio_all = ratio
        if dep:
            assert ratio >= 1/msqrt(2), (zz, ratio)
            if minratio_dep is None or ratio < minratio_dep[0]:
                minratio_dep = (ratio, H, zz)
        # in one unit class: even index => mod-2 deficient
        if I and I % 2 == 0:
            assert dep1
        if I % 2 == 1:
            assert not dep1
print("windows", nwin)
print("index distribution", sorted(idx.items()))
print("(dep_with_units, dep_exponents_only, one_unit_class) counts", sorted(par.items()))
print("min Delta*H^(1/4) over all windows", minratio_all)
print("min Delta*H^(1/4) over squareclass-dependent windows", minratio_dep[0] if minratio_dep else None)
