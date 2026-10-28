# Referee J4, check 2 (exact): J4 steps 1-2 on actual clusters, with the bookkeeping the judge's own check skipped.
# For a 5-subset S of lattice points of x^2+y^2 = N (N arbitrary: ramified 2, inert q, common split factors allowed):
#   g_S := gcd_{Z[i]} (Euclid), N'' := N/|g_S|^2.
#   (i)   N'' == prod_{p split} p^(max_S a(p) - min_S a(p))   with a(p) = v_{pi_p}
#   (ii)  formula gcd ~ (1+i)^{e0} prod q^{f_q} prod pi^{min a} pibar^{e_p - max a} (associate check)
#   (iii) v_{1+i}(z) = e0 and v_q(z) = f_q for every point (common ramified/inert part)
#   (iv)  zeta = z/g_S primitive (gcd unit), distinct, |zeta|^2 = N''; [zeta] = [z] (ratios equal)
#   (v)   |g_S| <= min chord |z_a - z_b|  (Lemma 1.2), so R'' >= R / chord
#   (vi)  the unit classes (eps_j) of the points in S are recorded: identity must not depend on them
import random, sys
from itertools import combinations
from math import atan2, sqrt, log
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from gz import *

def decompose(N):
    f = factor_int(N)
    e0 = f.get(2, 0)
    inert = {q: e for q, e in f.items() if q % 4 == 3}
    split = {p: e for p, e in f.items() if p % 4 == 1}
    return e0, inert, split

def unit_of(z, e0, inert, split, pis):
    w = z
    for _ in range(e0): w = gdivexact(w, (1, 1))
    for q, e in inert.items():
        for _ in range(e//2): w = gdivexact(w, (q, 0))
    for p, e in split.items():
        a = val(z, pis[p])
        for _ in range(a): w = gdivexact(w, pis[p])
        for _ in range(e-a): w = gdivexact(w, gconj(pis[p]))
    assert w is not None and norm(w) == 1, (z, w)
    return w

stats = dict(subsets=0, fail=0, cross_unit=0, ramified=0, inert=0, nonnormalised=0)
def check_subset(pts, N, e0, inert, split, pis, units):
    S = pts
    g = ggcd_list(S)
    assert N % norm(g) == 0
    N2 = N // norm(g)
    prod = 1
    gform = (1, 0)
    for _ in range(e0): gform = gmul(gform, (1, 1))
    for q, e in inert.items():
        for _ in range(e//2): gform = gmul(gform, (q, 0))
    for p, e in split.items():
        a = [val(z, pis[p]) for z in S]
        prod *= p ** (max(a) - min(a))
        for _ in range(min(a)): gform = gmul(gform, pis[p])
        for _ in range(e - max(a)): gform = gmul(gform, gconj(pis[p]))
    ok = (prod == N2)
    q1 = gdivexact(g, gform); ok &= (q1 is not None and is_unit(q1))
    for z in S:
        ok &= (val(z, (1, 1)) == e0)
        for q, e in inert.items(): ok &= (val(z, (q, 0)) == e//2)
    zeta = [gdivexact(z, g) for z in S]
    ok &= all(t is not None for t in zeta)
    ok &= is_unit(ggcd_list(zeta)) and len(set(zeta)) == 5 and all(norm(t) == N2 for t in zeta)
    # [zeta] = [z]: zeta_j * z_1 == z_j * zeta_1
    ok &= all(gmul(zeta[j], S[0]) == gmul(S[j], zeta[0]) for j in range(5))
    mind2 = min(norm((a[0]-b[0], a[1]-b[1])) for a, b in combinations(S, 2))
    ok &= norm(g) <= mind2
    stats['subsets'] += 1
    if not ok: stats['fail'] += 1; print("FAIL", N, S)
    if len(set(units[z] for z in S)) > 1: stats['cross_unit'] += 1
    if e0 > 0: stats['ramified'] += 1
    if inert: stats['inert'] += 1
    if norm(g) > 1: stats['nonnormalised'] += 1
    return N2

def run_circle(N, max_subsets=3000, rng=random.Random(1)):
    e0, inert, split = decompose(N)
    if any(e % 2 for e in inert.values()): return
    pis = {p: gauss_prime_above(p) for p in split}
    pts = lattice_points(N)
    if len(pts) < 5: return
    units = {z: unit_of(z, e0, inert, split, pis) for z in pts}
    allS = list(combinations(range(len(pts)), 5))
    if len(allS) > max_subsets: allS = rng.sample(allS, max_subsets)
    for S in allS: check_subset([pts[i] for i in S], N, e0, inert, split, pis, units)

if __name__ == "__main__":
    # (A) the record's 10-point cluster, and the same cluster times gamma = (1+i)^3 * 3 * 7 * (2+i)^2 (non-normalised)
    n = 1176852625
    pts = lattice_points(n)
    L = len(pts); best = None
    for i in range(L):
        cl = [pts[(i+j) % L] for j in range(10)]
        span = sqrt(norm((cl[0][0]-cl[-1][0], cl[0][1]-cl[-1][1])))
        if best is None or span < best[0]: best = (span, cl)
    cl = best[1]
    print("record cluster: n=%d, 10-point chord span/sqrt(R) = %.6f" % (n, best[0]/n**0.25))
    gamma = (1, 0)
    for f in [(1, 1)]*3 + [(3, 0), (7, 0), (2, 1), (2, 1)]: gamma = gmul(gamma, f)
    for tag, mult in (("normalised", (1, 0)), ("times gamma=(1+i)^3*3*7*(2+i)^2", gamma)):
        C2 = [gmul(z, mult) for z in cl]; N = norm(C2[0]); assert all(norm(z) == N for z in C2)
        e0, inert, split = decompose(N); pis = {p: gauss_prime_above(p) for p in split}
        units = {z: unit_of(z, e0, inert, split, pis) for z in C2}
        before = dict(stats)
        for S in combinations(C2, 5): check_subset(list(S), N, e0, inert, split, pis, units)
        print("  %s: N=%d, subsets=%d, fails=%d, unit classes in cluster=%d" %
              (tag, N, stats['subsets']-before['subsets'], stats['fail']-before['fail'], len(set(units.values()))))
    # (B) exhaustive over all N <= NMAX with >= 5 lattice points: all 5-subsets (sampled to 3000 per circle)
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    cnt = 0
    for N in range(1, NMAX+1):
        f = factor_int(N)
        r = 1
        for p, e in f.items():
            if p % 4 == 3 and e % 2: r = 0
            if p % 4 == 1: r *= e+1
        if 4*r < 5: continue
        run_circle(N, max_subsets=300); cnt += 1
    print("circles N<=%d with >=5 points: %d" % (NMAX, cnt))
    print("totals:", stats)
