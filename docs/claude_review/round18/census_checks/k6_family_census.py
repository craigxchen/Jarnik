"""k = 6 calibration (Data): the L8 Thm B family E_t at t0 = (2,3,-5,1), orbit Q1 + N(Q2 - Q1).
Uses the repository's k6fam scripts only to produce the rational six-point configurations
(x1, x2, x3, a, b, c); all contact data are recomputed with census_lib (factoring-free gcd method),
anchored at x1.  Reports w* , eps(Phi) (intrinsic), and eps after removing the structural bad set
Sigma (L8 §3.3), whose primes are moved to residues."""
import sys
sys.path.insert(0, 'k6fam')
from fractions import Fraction as Fr
from math import log
from lift import lift_T
from ell import quotient_curve
from census_lib import integral_chart, contact_integers_gcd, check_boolean_identity
SIGMA = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 61, 67, 71, 73, 97, 113, 311, 521, 571, 599, 1321, 2281]
xs = [Fr(2), Fr(3), Fr(-5), Fr(1)]
F, C, g, pts = quotient_curve(xs)
kk = g[3]
Q1, Q2 = pts[0], pts[1]
D = C.add(Q2, (Q1[0], -Q1[1]))
R = Q1
print("N   w*      eps     eps_SigmaFree  w_SigmaFree  log(Sigma-part total)  maxlog|X|")
for N in range(1, 15):
    R = C.add(R, D)
    L = lift_T(F, R[0] / kk)
    if not L:
        print(N, 'no lift'); continue
    pt = L[0]
    points = [xs[0], xs[1], xs[2], pt[0], pt[1], pt[2]]
    if len(set(points)) < 6:
        print(N, 'degenerate'); continue
    X = integral_chart(points, 0)
    d, delta = contact_integers_gcd(X)
    assert delta == 1
    lg = [log(v) if v > 1 else 0.0 for v in d.values()]
    lo, hi = min(lg), max(lg)
    eps = (hi - lo) / (hi + lo) if lo > 0 else 1.0
    # Sigma-free parts
    sf = {}
    sig = 0.0
    for T, v in d.items():
        u = v
        for p in SIGMA:
            while u % p == 0:
                u //= p
        sf[T] = u
        sig += log(v) - log(u)
    lg2 = [log(v) if v > 1 else 0.0 for v in sf.values()]
    lo2, hi2 = min(lg2), max(lg2)
    eps2 = (hi2 - lo2) / (hi2 + lo2) if lo2 > 0 else 1.0
    if N <= 8:
        assert check_boolean_identity(X) if max(abs(x) for x in X) < 10**60 else True
    print("%2d %8.2f  %.4f  %.4f        %8.2f      %7.2f              %.1f" % (
        N, (hi + lo) / 2, eps, eps2, (hi2 + lo2) / 2, sig, max(log(abs(x)) for x in X if x)))
