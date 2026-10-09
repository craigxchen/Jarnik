from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations
from lift import lift_T
from ell import quotient_curve
from profile import profile

xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
F, C, g, pts = quotient_curve(xs)
k = g[3]
Q1, Q2 = pts[0], pts[1]
D = C.add(Q2, (Q1[0], -Q1[1]))  # Q2 - Q1
R = Q1
bad = {}
for N in range(1, 15):
    R = C.add(R, D)          # Q1 + N (Q2 - Q1)
    L = lift_T(F, R[0] / k)
    if not L:
        print(N, 'no lift'); continue
    try:
        d, b, Dm = profile(xs, L[0])
    except ZeroDivisionError:
        print(N, 'boundary point'); continue
    logs = [log(x) for x in d.values()]
    mn, mx = min(logs), max(logs)
    shared = {}
    for A, B in combinations(d, 2):
        gg = gcd(d[A], d[B])
        if gg > 1:
            shared[gg] = shared.get(gg, 0) + 1
    print('N=%2d  w~mean=%9.2f  min/mean=%.3f  max/mean=%.3f  max b_ij=%d  shared gcds=%s' % (
        N, sum(logs) / 25, mn / (sum(logs) / 25), mx / (sum(logs) / 25), max(b.values()), sorted(shared)[:6]))
