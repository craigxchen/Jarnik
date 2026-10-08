# Referee: exhaustive exact checks over every circle x^2+y^2=N, N<=NMAX.
#  (i)  diameter claim: every 5 lattice points have diam >= (2^14*36)^(1/20) N^(1/5),
#       checked exactly as diam^20 >= 2^14*36*N^4 on all cyclic 5-windows (windows suffice:
#       a 5-set inside an open half-plane has diam = chord between its extreme points, which
#       is >= that of the 5-window starting at the first extreme point; a 5-set not inside an
#       open half-plane has diam >= sqrt(3) R which exceeds the bound for N>=2).
#  (ii) (B2): for the PRIMITIVE window, 128*sigma*prod s_ij * N^2 <= (dmax^2)^5 (exact),
#       sigma = product of p over singleton layers.
#  (iii) 6 | L_0 for every primitive 5-window inside an open half-plane; also random 5-subsets.
import sys, math, random
from fractions import Fraction
from math import gcd, isqrt
sys.path.insert(0, '.')
from rlib import gmul, gconj, gnorm, gdivmod_exact, ggcd_list, gauss_prime_above, vpi, angle

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
spf = list(range(NMAX + 1))
for i in range(2, isqrt(NMAX) + 1):
    if spf[i] == i:
        for j in range(i * i, NMAX + 1, i):
            if spf[j] == j: spf[j] = i
def factor(n):
    f = {}
    while n > 1:
        p = spf[n]; f[p] = f.get(p, 0) + 1; n //= p
    return f
pts = {}
for x in range(0, isqrt(NMAX) + 1):
    for y in range(0, isqrt(NMAX - x * x) + 1):
        n = x * x + y * y
        if n == 0: continue
        pts.setdefault(n, set()).update({(x, y), (-x, y), (x, -y), (-x, -y)})
BOUND = 2 ** 14 * 36
stats = dict(circles=0, windows=0, diam_ok=0, B2_ok=0, six_ok=0, halfplane=0, rnd=0, rnd_six=0)
min_ratio = None; tight_B2 = None
PI_CACHE = {}
def cluster_data(P):
    """P primitive list of 5 points in an open half plane. returns (N, sigma, prods, L0, dmax2)."""
    N = gnorm(P[0])
    prods = 1; L0 = 1; dmax2 = 0
    for i in range(5):
        for j in range(i + 1, 5):
            a, b = P[i], P[j]
            dot = a[0] * b[0] + a[1] * b[1]; det = a[0] * b[1] - a[1] * b[0]
            s = Fraction(N + dot, det).denominator
            prods *= s; L0 = L0 * s // gcd(L0, s)
            d2 = (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
            dmax2 = max(dmax2, d2)
    sigma = 1
    for p, e in factor(N).items():
        assert p % 4 == 1
        pi = PI_CACHE.setdefault(p, gauss_prime_above(p))
        lv = [vpi(z, pi) for z in P]
        for tau in range(1, e + 1):
            S = sum(1 for v in lv if v < tau)
            assert 1 <= S <= 4
            if S in (1, 4): sigma *= p
    return N, sigma, prods, L0, dmax2

def in_open_halfplane(Q):
    angs = sorted(angle(z) for z in Q)
    gaps = [(angs[(t + 1) % len(Q)] - angs[t]) % (2 * math.pi) for t in range(len(Q))]
    return max(gaps) > math.pi + 1e-9

random.seed(5)
for N in range(1, NMAX + 1):
    if N not in pts or len(pts[N]) < 5: continue
    stats['circles'] += 1
    L = sorted(pts[N], key=angle); n = len(L)
    for st in range(n):
        W = [L[(st + t) % n] for t in range(5)]
        stats['windows'] += 1
        d2 = (W[0][0] - W[4][0]) ** 2 + (W[0][1] - W[4][1]) ** 2
        if in_open_halfplane(W):
            stats['halfplane'] += 1
            dmax2 = max((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for a in W for b in W)
            if dmax2 ** 10 >= BOUND * N ** 4: stats['diam_ok'] += 1
            else: print('DIAMETER VIOLATION', N, W)
            r = (dmax2 ** 10 / N ** 4) ** (1 / 20)
            if min_ratio is None or r < min_ratio[0]: min_ratio = (r, N, W)
            g = ggcd_list(W); P = [gdivmod_exact(z, g) for z in W]
            Np, sigma, prods, L0, dm2 = cluster_data(P)
            if 128 * sigma * prods * Np ** 2 <= dm2 ** 5: stats['B2_ok'] += 1
            else: print('B2 VIOLATION', N, W)
            ratio = Fraction(dm2 ** 5, 128 * sigma * prods * Np ** 2)
            if tight_B2 is None or ratio < tight_B2[0]: tight_B2 = (ratio, N, W, sigma, prods)
            if L0 % 6 == 0: stats['six_ok'] += 1
            else: print('6|L0 VIOLATION', N, W, L0)
        else:
            # not in open half-plane: diameter >= sqrt(3)R; check claim directly anyway
            dmax2 = max((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for a in W for b in W)
            if dmax2 ** 10 >= BOUND * N ** 4: stats['diam_ok'] += 1
            else: print('DIAMETER VIOLATION (wide)', N, W)
    # random non-consecutive 5-subsets in a half-plane
    if n >= 8:
        for _ in range(3):
            Q = random.sample(L, 5)
            if not in_open_halfplane(Q): continue
            g = ggcd_list(Q); P = [gdivmod_exact(z, g) for z in Q]
            Np, sigma, prods, L0, dm2 = cluster_data(P)
            stats['rnd'] += 1
            if L0 % 6 == 0 and 128 * sigma * prods * Np ** 2 <= dm2 ** 5: stats['rnd_six'] += 1
            else: print('RANDOM SUBSET VIOLATION', N, Q, L0)
print(stats)
print('min diam/N^(1/5) over half-plane 5-windows: %.4f at N=%d %s' % min_ratio)
print('tightest (B2) ratio dmax^10/(128 sigma prod s N^2) = %.4f at N=%d %s sigma=%d prod s=%d'
      % (float(tight_B2[0]), tight_B2[1], tight_B2[2], tight_B2[3], tight_B2[4]))
