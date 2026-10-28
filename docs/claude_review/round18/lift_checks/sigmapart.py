from fractions import Fraction as Fr
from math import log
from ell import quotient_curve
from lift import lift_T
from anchored import anchored_form
SIGMA = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 61, 67, 71, 73, 97, 113, 311, 521, 571, 599, 1321, 2281]
def sigma_part(n):
    r = 1
    for p in SIGMA:
        while n % p == 0:
            n //= p; r *= p
    return r
xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
F, C, g, pts = quotient_curve(xs); k = g[3]
Q1, Q2 = pts[0], pts[1]
Dl = C.add(Q2, (Q1[0], -Q1[1])); R = Q1
for N in range(1, 15):
    R = C.add(R, Dl)
    L = lift_T(F, R[0] / k)
    try:
        Ymax, tmax, n = anchored_form(xs, L[0])
    except ZeroDivisionError:
        continue
    sp = [sigma_part(v) for v in n.values()]
    good = [v // s for v, s in zip(n.values(), sp)]
    lg = [log(x) for x in good]; w = sum(lg) / len(lg)
    print('N=%2d  total log(Sigma-parts) = %6.2f   good-prime blocks: min/w=%.3f max/w=%.3f  (w=%.1f)' % (N, sum(log(s) for s in sp), min(lg) / w, max(lg) / w, w))
