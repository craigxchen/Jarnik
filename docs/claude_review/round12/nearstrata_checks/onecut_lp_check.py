"""Consistency check of Theorem 2.5(a) with LP values (evidence): for one-cut boxes A_{S,L} (d <= 5), random u near strata,
eta in {0.3, 0.1, 0.03}: whenever nu_eta(u) >= m, the LP value kappa_V(u) = dist_inf(l_u^m, Poly_<m|_A)/d^m is
>= (m!/(2d)^m) eta^m. (kappa_V >= kappa_H, so this cannot fail if the theorem is right; it would expose a wrong proof.)"""
import random, math
from fractions import Fraction as F
from duality_check import box, primal
from one_cut_gadgets import build

rng = random.Random(7)
cnt = tested = 0; minratio = float('inf')
for trial in range(400):
    d = rng.choice([4, 5]); s = rng.randint(2, d); L = rng.randint(1, s - 1)
    S = list(range(s)); R = list(range(s, d)); n = d + 1
    m = rng.randint(2, 6)
    A = box(d, [(S, L)])
    base = [rng.choice([0.0, 0.0, 0.5, -0.4, rng.uniform(-1, 1)]) for _ in range(d)]
    u = [b + rng.uniform(-0.05, 0.05) for b in base]
    mx = max(abs(x) for x in u); u = [x / mx for x in u]
    for eta in (0.3, 0.1, 0.03):
        uF = [F(x).limit_denominator(10**9) for x in u]
        G, info = build(uF, S, R, L, F(eta).limit_denominator(1000))
        if info['nu'] >= m:
            tested += 1
            pv, _, _ = primal(A, m, u)
            kappa = pv / d**m
            bound = math.factorial(m) / (2 * d)**m * eta**m
            minratio = min(minratio, kappa / bound)
            assert kappa >= bound * (1 - 1e-7), (kappa, bound, d, s, L, m, u, eta)
print("tested", tested, "cases with nu_eta >= m; min kappa_V / bound =", minratio)
