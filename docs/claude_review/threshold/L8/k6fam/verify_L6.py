"""One-shot exact verification for Theorem B ((L_6) is false).  Run from this directory: python3 verify_L6.py
(1) open conditions defining U, at t0=(2,3,-5,1) and at random rational points;
(2) Q_1 non-torsion on E'_t (Mazur: nQ != O for n in 1..10,12) at the same points;
(3) dominance: rank of d(B_C)/dt = 2 for C={x_j,a};
(4) actual configurations on E_t0 along Q1 + N(Q2-Q1): balance ratio, Y_j = +-1, Sigma-supported overlaps."""
import random
from fractions import Fraction as Fr
from math import gcd, log
from itertools import combinations
from anchored import U_conditions, anchored_form
from ell import quotient_curve, nontorsion
from dominance import Du, contact_coords, rank
from lift import lift_T

def check_point(xs):
    conds, s, t = U_conditions(xs)
    inU = all(conds.values())
    F, C, g, pts = quotient_curve(xs)
    nt = inU and pts[0] is not None and nontorsion(C, pts[0])
    x = [Du(xs[0], [1, 0, 0]), Du(xs[1], [0, 1, 0]), Du(xs[2], [0, 0, 1]), Du(xs[3], [0, 0, 0])]
    rk = min(rank([u.g, v.g]) for u, v in (contact_coords(x, j) for j in range(3))) if inU else None
    return inU, nt, rk

print('(1)-(3) at t0 = (2,3,-5,1):', check_point([Fr(2), Fr(3), Fr(-5), Fr(1)]))
random.seed(2026)
stats = [0, 0, 0, 0]
for _ in range(40):
    xs = [Fr(random.randint(-60, 60), random.randint(1, 9)) for _ in range(3)] + [Fr(1)]
    try:
        inU, nt, rk = check_point(xs)
    except ZeroDivisionError:
        inU, nt, rk = False, False, None
    stats[0] += 1; stats[1] += inU; stats[2] += bool(nt); stats[3] += (rk == 2)
print('(1)-(3) random rational t (x4=1): tried %d, in U %d, Q1 non-torsion %d, dominance rank 2 %d' % tuple(stats))

xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
F, C, g, pts = quotient_curve(xs); k = g[3]
Q1, Q2 = pts[0], pts[1]
Dl = C.add(Q2, (Q1[0], -Q1[1])); R = Q1
SIGMA = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 61, 67, 71, 73, 97, 113, 311, 521, 571, 599, 1321, 2281}
def sig_part(n):
    r = 1
    for p in SIGMA:
        while n % p == 0: n //= p; r *= p
    return r
print('(4) orbit on E_t0:')
for N in range(1, 13):
    R = C.add(R, Dl)
    L = lift_T(F, R[0] / k)
    try:
        Ymax, tmax, n = anchored_form(xs, L[0])
    except ZeroDivisionError:
        continue
    good = {T: v // sig_part(v) for T, v in n.items()}
    lg = [log(v) for v in good.values()]; w = sum(lg) / len(lg)
    overlap = max(gcd(a, b) for a, b in combinations(good.values(), 2))
    print('  N=%2d  w=%7.1f  min/w=%.3f max/w=%.3f  max|Y_j|=%d  max gcd of Sigma-free blocks=%d  log Sigma-part total=%.1f' % (
        N, w, min(lg) / w, max(lg) / w, Ymax, overlap, sum(log(sig_part(v)) for v in n.values())))
