"""Brute force: for every N <= Nmax that is a sum of two squares, the minimum over
5-subsets of distinct lattice points on x^2+y^2=N of diam/N^(1/5), compared with
the bound 1.943 of Theorem B (five points: diam >= (2^14*36)^(1/20) R^(2/5))."""
import math, sys
from itertools import combinations
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
best = (1e9, None)
bound = (2 ** 14 * 36) ** (1 / 20)
for N in range(1, Nmax + 1):
    pts = []
    r = math.isqrt(N)
    for x in range(-r, r + 1):
        y2 = N - x * x
        y = math.isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    if len(pts) < 5:
        continue
    pts.sort(key=lambda p: math.atan2(p[1], p[0]))
    n = len(pts)
    # minimal-diameter 5-subsets are 5 angularly consecutive points
    for s in range(n):
        cl = [pts[(s + t) % n] for t in range(5)]
        d = max(math.dist(a, b) for a, b in combinations(cl, 2))
        k = d / N ** 0.2
        if k < best[0]:
            best = (k, (N, cl))
print('bound %.4f ; min diam/N^(1/5) over N<=%d: %.4f at %s' % (bound, Nmax, best[0], best[1]))
assert best[0] >= bound
