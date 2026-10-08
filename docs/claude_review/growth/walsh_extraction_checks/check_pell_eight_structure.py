"""Item 511 eight-point translated Pell family: verify exactly that the sign words form an affine 3-cube
(pure Sylvester core of order 8 on the four heavy blocks) and that the light F=1+2i class has the
non-affine majority pattern; report the normalised arc constant with exact chords."""
import math, itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gauss_common import *
from math import isqrt
def lam_pow(r):
    x, y = 1, 0
    a, b = (9, 4) if r >= 0 else (9, -4)
    for _ in range(abs(r)):
        x, y = x * a + 5 * y * b, x * b + y * a
    return x, y
F = (1, 2)
def H(r):
    x, y = lam_pow(r)
    return (x + y * F[0], y * F[1])
for n in (1, 61):
    blocks = [H(r) for r in (n - 2, n, n + 2, n + 4)]
    pts = []; words = []
    for w in itertools.product((0, 1), repeat=4):
        if sum(w) % 2 == 1:
            z = F if sum(w) == 1 else gconj(F)
            for bit, B in zip(w, blocks):
                z = gmul(z, B if bit else gconj(B))
            pts.append(z); words.append(w)
    N = gnorm(pts[0]); assert all(gnorm(z) == N for z in pts) and len(set(pts)) == 8
    g = pts[0]
    for z in pts[1:]: g = ggcd(g, z)
    P = [gdivexact(z, g) for z in pts]
    N2 = gnorm(P[0])
    # rotate all points by units into one quadrant-arc and take the largest exact chord^2
    best = None
    for us in itertools.product(range(4), repeat=7):
        Q = [P[0]] + [gmul(UNITS[u], z) for u, z in zip(us, P[1:])]
        mx = max(gnorm((a[0] - b[0], a[1] - b[1])) for a in Q for b in Q)
        if best is None or mx < best: best = mx
    # C >= chord/sqrt(R) = sqrt(max chord^2) / N2^(1/4)
    logC = 0.5 * math.log(best) - math.log(N2) / 4
    X = [int("".join(map(str, w)), 2) for w in words]
    D = set(x ^ X[0] for x in X); assert all((u ^ v) in D for u in D for v in D)
    maj = [1 if sum(w) == 1 else 0 for w in words]
    coords = [w[:3] for w in words]
    aff = any(all((sum(ai * wi for ai, wi in zip(a, cw)) + c) % 2 == m for cw, m in zip(coords, maj))
              for a in itertools.product((0, 1), repeat=3) for c in (0, 1))
    print("n=%d: 8 distinct equal-norm primitive points, log10 N=%d, log10(max chord/sqrt R)=%.2f ; block signs = affine 3-cube ; F-class affine: %s"
          % (n, len(str(N2)) - 1, logC / math.log(10), aff))
