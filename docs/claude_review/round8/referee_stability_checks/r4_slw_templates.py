#!/usr/bin/env python3
"""R4. Is the hypothesis behind stability.md Section 8.6's pointer "by weakest.md (draft), five simultaneous
forms do" consistent?  weakest.md's SLW(k, g, kappa) (with g < 1/4, k = 5) is refuted by vertical templates
(an independent re-derivation of the construction in referee_weakest_checks/slw_template_counterexample.py,
no factoring needed):

  b_1..b_D distinct positive integers, sigma_1..sigma_k in Z^D independent with
  sum_l sigma_il b_l^(2j-1) = 0 (j = 1..s), k = D - s.  For t divisible by L = lcm(b) * P^2
  (P = product of primes <= 2 max b), u_l = t/b_l + i are pairwise coprime to each other and to the
  conjugates, and beta_i = prod_l u_l^(sigma_il) (conj for negative exponents) has
     Lambda_i = Log(beta_i / conj beta_i) = 2 i sum_l sigma_il arctan(b_l / t) = O(t^-(2s+1)),
  while prod_{p in U} H(alpha_p) <= 2^n prod_l N(u_l) ~ t^(2D).
So gamma_N := log(1/max|Lambda_i|) / sum_l log N(u_l) -> (2s+1)/(2D) = (2s+1)/(2(k+s)), which exceeds 1/4
iff s > (k-2)/3 (k = 5: s >= 2), and tends to 1 as s -> infinity for fixed k.  Since n = sum omega(N(u_l)) = o(log t), the
factors c^(n+1) prod (1+B_j)^(-kappa) 2^(-g n) are t^(-o(1)); hence SLW(k, g, kappa) fails for every g < gamma.
Exact: the moment conditions (integers) and pairwise coprimality of the u_l and conj u_m (Gaussian gcd).
Angles: mpmath at 2000 digits.  This is evidence for the asymptotic statement plus an exact construction.
"""
from math import gcd
from fractions import Fraction
import mpmath
import sympy
from rlib import ggcd, gnorm, gconj

mpmath.mp.dps = 2000


def kernel_int(b, s):
    Mx = sympy.Matrix([[bl ** (2 * j - 1) for bl in b] for j in range(1, s + 1)])
    vecs = []
    for v in Mx.nullspace():
        den = sympy.ilcm(*[sympy.fraction(x)[1] for x in v])
        w = [int(x * den) for x in v]
        g = 0
        for x in w:
            g = gcd(g, x)
        vecs.append([x // g for x in w])
    return vecs


def run(b, s, Ts):
    D = len(b)
    sig = kernel_int(b, s)
    k = len(sig)
    for v in sig:
        for j in range(1, s + 1):
            assert sum(v[l] * b[l] ** (2 * j - 1) for l in range(D)) == 0
    assert sympy.Matrix(sig).rank() == k
    L = 1
    for x in b:
        L = L * x // gcd(L, x)
    P = 1
    for p in sympy.primerange(2, 2 * max(b) + 1):
        P *= p
    Q = L * P * P
    pred = Fraction(2 * s + 1, 2 * D)
    print(f"b = {b}, s = {s}: k = {k} forms, D = {D}, max|sigma| = {max(abs(x) for v in sig for x in v)}, "
          f"predicted gamma -> (2s+1)/(2D) = {pred} = {float(pred):.4f}")
    for T in Ts:
        t = Q * T
        us = [(t // bl, 1) for bl in b]
        cop = all(gnorm(ggcd(us[a], us[c])) == 1 and gnorm(ggcd(us[a], gconj(us[c]))) == 1
                  for a in range(D) for c in range(a))
        lam = [abs(2 * sum(v[l] * mpmath.atan(mpmath.mpf(b[l]) / t) for l in range(D))) for v in sig]
        mx = max(lam)
        hsum = sum(mpmath.log(gnorm(u)) for u in us)
        gam = -mpmath.log(mx) / hsum
        print(f"   T = {T:>14d}: log t = {float(mpmath.log(t)):8.2f}, coprime u's: {cop}, "
              f"max|Lambda| = {mpmath.nstr(mx, 5)}, gamma_N = {mpmath.nstr(gam, 6)}")
    return pred


def main():
    Ts = [10 ** 3 + 9, 10 ** 6 + 3, 10 ** 12 + 39, 10 ** 24 + 7, 10 ** 48 + 75]
    p1 = run([1, 2, 3, 4, 5, 6, 7, 8, 9], 4, Ts)     # k = 5 forms, gamma -> 1/2
    p2 = run(list(range(1, 12)), 6, Ts)             # k = 5 forms, gamma -> 13/22
    print("Five forms (k = 5) with gamma_N -> 1/2 and 13/22, both > 1/4: SLW(5, g, kappa) with g < 1/4 is "
          "violated for large t (any c, kappa), as are all SLW(5, g) with g < 13/22.")
    print("ALL CHECKS PASSED:", p1 > Fraction(1, 4) and p2 > Fraction(1, 4))


if __name__ == "__main__":
    main()
