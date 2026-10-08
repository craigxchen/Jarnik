# Referee cross-check of the author's CG-regime data: list every primitive 5-window (all N<=NMAX)
# with diam/N^(1/5) < KMAX, together with its span (A>=L iff span<=90 deg) and intrinsic L_0.
import sys, math
from fractions import Fraction
from math import gcd, isqrt
sys.path.insert(0, '.')
from rlib import gnorm, gdivmod_exact, ggcd_list, angle

NMAX = int(sys.argv[1]); KMAX = float(sys.argv[2])
pts = {}
for x in range(0, isqrt(NMAX) + 1):
    for y in range(0, isqrt(NMAX - x * x) + 1):
        n = x * x + y * y
        if n: pts.setdefault(n, set()).update({(x, y), (-x, y), (x, -y), (-x, -y)})
found = {}
for N in sorted(pts):
    L = sorted(pts[N], key=angle); n = len(L)
    if n < 5: continue
    for st in range(n):
        W = [L[(st + t) % n] for t in range(5)]
        span = (angle(W[4]) - angle(W[0])) % (2 * math.pi)
        if span >= math.pi - 1e-9: continue
        g = ggcd_list(W); P = [gdivmod_exact(z, g) for z in W]
        Np = gnorm(P[0])
        d = max(math.dist(a, b) for a in P for b in P)
        r = d / Np ** 0.2
        if r >= KMAX: continue
        L0 = 1
        for i in range(5):
            for j in range(i + 1, 5):
                a, b = P[i], P[j]
                s = Fraction(Np + a[0] * b[0] + a[1] * b[1], a[0] * b[1] - a[1] * b[0]).denominator
                L0 = L0 * s // gcd(L0, s)
        key = (Np, round(r, 6))
        found[key] = (r, Np, L0, math.degrees(span))
rows = sorted(found.values())
print('primitive 5-windows with diam/N^(1/5) <', KMAX, 'and N<=', NMAX, ':', len(rows))
for r, Np, L0, sp in rows[:25]:
    print('  K5=%.4f N=%d L0=%d span=%.1f deg %s' % (r, Np, L0, sp, 'A>=L' if sp <= 90 + 1e-9 else 'A<L'))
