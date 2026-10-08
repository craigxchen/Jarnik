"""Exact discriminant squareclasses on both isotropic rulings."""

from fractions import Fraction as F
from itertools import combinations
from math import isqrt, lcm, prod

from check_six_point_isotropic_circle_cover import determinant, probe, pscale
from check_segre_gradient_arithmetic import rref


def rational_square(value):
    value = F(value)
    return (value > 0 and isqrt(value.numerator) ** 2 == value.numerator
            and isqrt(value.denominator) ** 2 == value.denominator)


def dual_frame(u, rows):
    dot = lambda a, b: sum(c*x*y for c, x, y in zip(u, a, b))
    for inds in combinations(range(6), 3):
        matrix = [[rows[i][j]*u[j] for j in inds]
                  + [F(i == k) for k in range(3)] for i in range(3)]
        reduced, pivots = rref(matrix, 3)
        if len(pivots) == 3:
            break
    temporary = [[F(0)]*6 for _ in range(3)]
    for i in range(3):
        for k, j in enumerate(inds):
            temporary[i][j] = reduced[k][3+i]
    gram = [[dot(a, b) for b in temporary] for a in temporary]
    return [[temporary[i][j]-sum(gram[i][k]*rows[k][j]/2
                                for k in range(3))
             for j in range(6)] for i in range(3)]


def primes(limit):
    return [p for p in range(3, limit)
            if all(p % d for d in range(2, isqrt(p)+1))]


def check(u, rows, dual):
    h = lcm(*(F(c).denominator for row in rows+dual for c in row))
    counts = [0, 0, 0]  # quadratics, rationally split quadratics, prime checks
    for ruling in (0, 1):
        if ruling == 0:
            x = [[rows[1][j], dual[2][j]] for j in range(6)]
            y = [[rows[2][j], -dual[1][j]] for j in range(6)]
        else:
            x = [[rows[1][j], rows[2][j]] for j in range(6)]
            y = [[dual[2][j], -dual[1][j]] for j in range(6)]
        x, y = ([pscale(v, h) for v in family] for family in (x, y))
        for inds in combinations(range(6), 3):
            subtotal = sum(u[j] for j in inds)
            if subtotal == 0:
                continue
            target = -subtotal*prod(u[j] for j in inds)
            coefficients = determinant([[[1], x[j], y[j]] for j in inds])
            coefficients += [F(0)]*(3-len(coefficients))
            assert all(F(c).denominator == 1 for c in coefficients)
            a, b, c = map(int, coefficients)
            discriminant = b*b-4*a*c
            assert discriminant != 0
            assert rational_square(F(discriminant, target))
            counts[0] += 1
            counts[1] += int(rational_square(target))
            exceptional = 2*h*subtotal*prod(u)
            for p in primes(180):
                if exceptional % p == 0:
                    continue
                assert discriminant % p != 0
                roots = sum((a+b*r+c*r*r) % p == 0 for r in range(p))
                roots += int(c % p == 0)  # parameter infinity
                legendre = pow(target % p, (p-1)//2, p)
                assert roots == (2 if legendre == 1 else 0)
                counts[2] += 1
    return counts


if __name__ == '__main__':
    totals = [0, 0, 0]
    for parameters in ([0, 1, 2, 3, 4, 5], [0, 1, 2, 4, 7, 11]):
        u, rows, dual, _, _ = probe(parameters)
        totals = [a+b for a, b in zip(totals, check(u, rows, dual))]
    # Resonance elsewhere is allowed by this selected-partition theorem.
    # This example includes rationally split quadratic factors and infinity roots.
    u = [1, -1, 1, -1, 1, -1]
    rows = [[F(1)]*6, list(map(F, [1, 1, 0, 0, 0, 0])),
            list(map(F, [0, 0, 1, 1, 0, 0]))]
    totals = [a+b for a, b in zip(totals, check(u, rows, dual_frame(u, rows)))]
    assert totals[1] > 0
    print('PASS:', totals[0], 'quadratics on both rulings;', totals[1],
          'rationally split cases;', totals[2], 'good-prime splitting checks.')
