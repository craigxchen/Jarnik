"""Corollaries of (IV.2) for the record's model M(D,c) (exact rational arithmetic).
Dictionary: L = 2D, E_ij = 0, delta = floor(D/2+c) - D/2 (any sign allowed: (IV.3) is a plain sum).
An m-point cluster (any m-subset of a cluster is a cluster) must satisfy
   (3/4)(m-24)(2D) + Sigma <= m(m+3) delta ,  Sigma >= 0.
"""
from fractions import Fraction as Fr
import math, sympy as sp

# 1. the final substitution of Part IV, symbolically
s, m = sp.symbols('s m')
lhs = sp.Rational(2, 9) * (m + 3) * (s - m / 4) - sp.Rational(2, 9) * s * (s - 1) + 4 - sp.Rational(8, 9) * s
rhs = 4 - m / 6 - sp.Rational(2, 9) * (s - m / 2) ** 2
assert sp.simplify(lhs - rhs) == 0
# coefficient of S: 8/9 + (4/9)(m-1)/2 = 2(m+3)/9
assert sp.simplify(sp.Rational(8, 9) + sp.Rational(4, 9) * (m - 1) / 2 - 2 * (m + 3) / 9) == 0
print("per-layer identity and S-coefficient: OK")


def delta(D, c):
    return Fr(math.floor(Fr(D, 2) + c)) - Fr(D, 2)


def max_m(D, c):
    """bound on M(D,c) from (IV.2) applied to m-subsets: (smallest excluded m) - 1; None if none excluded."""
    dl = delta(D, c)
    for mm in range(3, 3000):
        if Fr(3, 4) * (mm - 24) * 2 * D > mm * (mm + 3) * dl:
            return mm - 1
    return None


# 2. c = 0
for D in range(1, 120):
    mm = max_m(D, 0)
    assert mm <= 24
    if D % 2 == 1:
        assert mm <= 23
print("c=0: M(D,0) <= 24 for D<120 (all D by the formula); odd D: <= 23;",
      "odd-D threshold for 23 points: D >=", Fr(23 * 26, 3), "; for 22:", Fr(22 * 25, 6))

# 3. thresholds for c > 0: m-point clusters (m >= 25) impossible if D > f(m) c, f = 2m(m+3)/(3(m-24))
f = {mm: Fr(2 * mm * (mm + 3), 3 * (mm - 24)) for mm in range(25, 200)}
mbest = min(f, key=lambda k: f[k])
print("min f(m) at m =", mbest, "f =", f[mbest], float(f[mbest]), "; f(48) =", f[48], "; f(25) =", f[25], float(f[25]))
# so: D > f(m) c  =>  M(D,c) <= m-1
for c in (Fr(1, 2), 1, 2, 3):
    rows = []
    for D in (10, 34, 35, 68, 69, 136, 137, 467, 468, 934, 935, 1000):
        rows.append((D, max_m(D, c) if max_m(D, c) is not None else None))
    print("c =", c, [(D, mm if mm is not None else 'none') for D, mm in rows])
# 4. uniform bound sup_D M(D,c) <= max over D <= 68c of record bound 2^(2c+2)(D+2), and <= 47 beyond.
for c in (1, 2, 3):
    print("c=%d: sup_D M(D,c) <= max(47, %d)" % (c, 2 ** (2 * c + 2) * (68 * c + 2)))
