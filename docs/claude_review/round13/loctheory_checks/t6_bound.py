"""Rigorous all-j bound for t = 6 (b = 63): rho(Y_6(j)) > 9j + 9 for all j >= j0.
N(e) <= binom(e+b, b);  |Y_6(j)| >= (2j-30)^63 / 2^129  (Prop 5.5 with floor((2j+1)/32) >= (2j-30)/32).
The ratio (9j+72)/(2j-30) decreases in j, so  binom(9j+72,63) <= (9j+72)^63/63! < (2j-30)^63/2^129  at j0 implies all j >= j0.
Checked with exact integers at j0, plus the exact (non-relaxed) comparison at j0."""
from math import comb, factorial
b = 63
def ok_relaxed(j):
    return (9 * j + 72) ** 63 * 2 ** 129 < (2 * j - 30) ** 63 * factorial(63)
j0 = next(j for j in range(16, 10 ** 6) if ok_relaxed(j))
print('first j with relaxed inequality:', j0, ' holds at j0, j0+1, 10^4:', ok_relaxed(j0), ok_relaxed(j0 + 1), ok_relaxed(10 ** 4))
def exact_ok(j):
    Ylow = (2 ** (5 * b) // 2 ** 129) * ((2 * j + 1) // 32) ** b
    return comb(9 * j + 9 + b, b) < Ylow
print('exact (unrelaxed) check at j0..j0+5:', [exact_ok(j) for j in range(j0, j0 + 6)])
