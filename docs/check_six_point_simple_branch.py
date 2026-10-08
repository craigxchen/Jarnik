"""Exact first-order conic check at rational parabolic branches."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, lcm, prod
from random import Random

from check_segre_gradient_arithmetic import rref


def primitive_weights(ts):
    raw = [F(1, prod(t - s for j, s in enumerate(ts) if j != i))
           for i, t in enumerate(ts)]
    denominator = lcm(*(c.denominator for c in raw))
    values = [int(c * denominator) for c in raw]
    common = reduce(gcd, values)
    return [c // common for c in values]


def interpolate(xs, ys):
    matrix = [[F(t ** j) for j in range(5)] + [F(y)] for t, y in zip(xs, ys)]
    rows, pivots = rref(matrix, 5)
    assert len(pivots) == 5
    coefficients = [rows[j][5] for j in range(5)]
    assert all(sum(c * t ** j for j, c in enumerate(coefficients)) == y
               for t, y in zip(xs, ys))
    return coefficients


def check(ts):
    u = primitive_weights(ts)
    e = [[F(1)] * 6, list(map(F, ts)), [F(t * t) for t in ts]]
    beta = lambda a, b: sum(c * x * y for c, x, y in zip(u, a, b))
    assert all(beta(a, b) == 0 for a in e for b in e)
    # Any three columns of (1,t,t^2) are independent.
    matrix = [[e[i][j] * u[j] for j in range(3)] + [F(i == k) for k in range(3)]
              for i in range(3)]
    reduced, pivots = rref(matrix, 3)
    assert len(pivots) == 3
    temp = [[reduced[j][3 + k] if j < 3 else F(0) for j in range(6)] for k in range(3)]
    gram = [[beta(a, b) for b in temp] for a in temp]
    dual = [[temp[i][j] - sum(gram[i][k] * e[k][j] / 2 for k in range(3))
             for j in range(6)] for i in range(3)]
    assert all(beta(e[i], dual[j]) == int(i == j) for i in range(3) for j in range(3))
    assert all(beta(a, b) == 0 for a in dual for b in dual)
    dx, dy = dual[2], [-c for c in dual[1]]
    rhs = [2 * t * a - b for t, a, b in zip(ts, dx, dy)]
    derivative = interpolate(ts, rhs)
    s5 = sum(c * t ** 5 for c, t in zip(u, ts))
    assert s5 and derivative[4] * s5 == 3
    assert -4 * derivative[4] == F(-12, s5)
    no_subsum = all(sum(u[i] for i in inds) != 0
                    for size in (1, 2, 3) for inds in combinations(range(6), size))
    return no_subsum


def main():
    rng = Random(20260913)
    samples = [list(range(6)), [0, 1, 2, 4, 7, 11]]
    samples.extend(sorted(rng.sample(range(-20, 31), 6)) for _ in range(100))
    counts = [0, 0]
    for ts in samples:
        counts[check(ts)] += 1
    assert all(counts)
    print(f'PASS: {sum(counts)} exact branch derivatives; {counts[1]} weights without proper vanishing subsums.')
    print("Verified C' S5=3 and Delta'=-12/S5 in all cases.")


if __name__ == '__main__':
    main()
