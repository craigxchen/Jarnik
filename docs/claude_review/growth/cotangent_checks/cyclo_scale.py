"""Intrinsic cotangent scale L_0 of the five-row cyclotomic family.

Rebuilds the five Gaussian rows of cyclotomic_five_row_subendpoint_family.md
(k odd >= 19, k = 0,1 mod 3; d odd) and computes, without Gaussian gcds,
  c_ij = (N + <Z_i,Z_j>)/det(Z_i,Z_j)   (scale invariant),
  L_0  = lcm of reduced denominators,
  log N_prim via the all-edge lcm of reduced edge norms.
"""
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd, log


def fib_pair(n):
    if n == 0:
        return 0, 1
    a, b = fib_pair(n // 2)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def gfib(m):
    r = (m - 1) // 2
    a, b = fib_pair(r)
    return (b, a)


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def conj(a):
    return (a[0], -a[1])


def rows(k, d):
    x, y, z, w = k - 2, k, k + 2, k + 4
    ns = [x * z * w, y * y * w, y * z * z, y * z * w]
    H = [gfib(n * d) for n in ns]
    F, Fb = (1, 2), (1, -2)
    pats = ((0,), (1,), (2,), (3,), (0, 1, 2))
    Z = []
    for j, fw in enumerate(pats):
        v = Fb if j == 4 else F
        for t in range(4):
            v = mul(v, H[t] if t in fw else conj(H[t]))
        Z.append(v)
    return Z


def lcm(a, b):
    return a // gcd(a, b) * b


def analyse(k, d):
    Z = rows(k, d)
    N = Z[0][0] ** 2 + Z[0][1] ** 2
    cots = {}
    for i, j in combinations(range(5), 2):
        dot = Z[i][0] * Z[j][0] + Z[i][1] * Z[j][1]
        det = Z[i][0] * Z[j][1] - Z[i][1] * Z[j][0]
        cots[(i, j)] = Fraction(N + dot, det)
    L0 = 1
    for c in cots.values():
        L0 = lcm(L0, c.denominator)
    # primitive squared radius = lcm of reduced edge norms (a^2+s^2)/eps
    Np = 1
    for c in cots.values():
        a, s = c.numerator, c.denominator
        eps = 2 if (a % 2 and s % 2) else 1
        Np = lcm(Np, (a * a + s * s) // eps)
    # largest chord^2 = 4 Np s^2/(eps n) over pairs (primitive circle)
    best = 0.0
    lg = log(Np)
    for c in cots.values():
        a, s = c.numerator, c.denominator
        eps = 2 if (a % 2 and s % 2) else 1
        n = (a * a + s * s) // eps
        lchord2 = log(4) + lg + 2 * log(s) - log(eps) - log(n)
        best = max(best, lchord2)
    # log of max chord / R^(1/2) and / R^(2/5)
    logCstar = best / 2 - lg / 4
    logK5 = best / 2 - lg / 5
    ss = sorted(c.denominator for c in cots.values())
    print('k=%d d=%d: log N_prim=%.1f  log L0=%.2f  log C*=%.3f  '
          'log(maxchord/R^(2/5))=%.1f  residues=%s' %
          (k, d, lg, log(L0), logCstar, logK5, [round(log(s),1) for s in ss]))


if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 19
    for d in (1, 3, 5, 7, 9, 11):
        analyse(k, d)
