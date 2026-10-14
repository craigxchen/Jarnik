"""Corollary 3.7 of nearstrata.md: separation sigma_F (Cor 3.5) versus the smallest radius rho_F that a near-strata
statement can have near a cheap span, given the proved transversal order o_F >= m' - d_F (Lemma 5.0):
   kappa_H <= 4^{m'} eta^{o_F} at distance eta  =>  any valid radius satisfies 4^{m'} rho_F^{o_F} >= eps_0, i.e.
   rho_F >= rho_F^min := (eps_0 4^{-m'})^{1/o_F}.
sigma_F > rho_F^min for every cheap F iff for all (s, e):
   o (Q abar/2 + ln 2) < W/4 + m' ln 4 - ln(C d/(m'+1)),   o = m' - 2js - 2e (worst case of the proved order),
   Q = R(s,e) = 2^{s+1} N_2/((n-e)(n-e-2^s)),  abar = W/(2d) + 2 ln C'.
Linear in W: holds for all W >= W0 iff it holds at W0 and the W-coefficient 1/4 - o Q/(4d) is > 0 (condition B of compat.py).
We take C = 1, W0 = W_min(n) = 3 n ln n - 14 n (two_thirds (T3), K_0 < 13.91), and report the last failing n (step 1 to 2*10^4)."""
import math, numpy as np
from compat import params

def margin(n, r, W, C=1.0):
    d, w, j, mp = params(n, r)
    if w < 10 or 2 * j > n: return None
    N2 = n * (n - 1) / 2; Cp = max(C, 1.0)
    abar = W / (2 * d) + 2 * math.log(Cp)
    b = (mp - 1) // (2 * j); worst = float('inf'); worst_coef = float('inf')
    for s in range(b + 1):
        emax = (mp - 1) // 2 - j * s
        if emax < 0: continue
        e = np.arange(emax + 1, dtype=float); q = 2.0 ** s
        Q = 2 * q * N2 / ((n - e) * (n - e - q))
        o = np.maximum(mp - 2 * j * s - 2 * e, 0)
        G = W / 4 + mp * math.log(4) - math.log(C * d / (mp + 1)) - o * (Q * abar / 2 + math.log(2))
        worst = min(worst, float(G.min()))
        worst_coef = min(worst_coef, float((0.25 - o * Q / (4 * d)).min()))
    return worst, worst_coef

for name, rf in [('r = n', lambda n: n), ('r = n ln n', lambda n: n * math.log(n)), ('r = n^2', lambda n: n * n)]:
    last = None
    for n in range(500, 20001):
        Wmin = 3 * n * math.log(n) - 14 * n
        res = margin(n, rf(n), Wmin)
        if res is None or res[0] <= 0 or res[1] <= 0: last = n
    print(f"{name}: last n in [500, 20000] where sigma_F <= rho_F^min for some cheap F at some W >= W_min: {last}")
for n in (3000, 5185, 10**4):
    print(n, [round(x, 3) for x in margin(n, n, 3 * n * math.log(n) - 14 * n)])
