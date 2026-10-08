"""Exact checks for the fixed six-point moment elliptic obstruction."""

from fractions import Fraction as F
from itertools import combinations
from math import isqrt, lcm

from check_integer_cotangent_spectral import (
    ONE, add, conj, div, mul, neg, power, primitive_vector,
)
from check_integer_cotangent_normalization import edge_norm


P = (F(1), F(14))
U = (5, -5, 5, -5, -1, 1)


def ec_add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    a, b = q
    if x == a:
        if y == -b:
            return None
        slope = (3*x*x-204*x+297)/(2*y)
    else:
        slope = (b-y)/(a-x)
    c = slope*slope+102-x-a
    return c, -y+slope*(x-c)


def rational_sqrt(x):
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    assert a*a == x.numerator and b*b == x.denominator
    return F(a, b)


def ordinary_valuation(value, prime):
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def gaussian_valuation(z, prime):
    exponent = 0
    while True:
        quotient = div(z, prime)
        if any(F(c).denominator != 1 for c in quotient):
            return exponent
        z = tuple(int(c) for c in quotient)
        exponent += 1


def configuration(point):
    x, y = point
    assert y*y == x*(x-3)*(x-99)
    q = rational_sqrt(x)
    p = y/(99-x)
    assert p*p*(99-q*q) == q*q*(3-q*q)
    first = (1-p*q, p+q)
    second = (1+p*q, p-q)
    third = (1-q*q, 10*p)
    rows = [first, conj(first), second, conj(second), third, conj(third)]
    assert len(set(rows)) == 6
    norms = [a*a+b*b for a, b in rows]
    assert len(set(norms)) == 1
    for exponent in range(3):
        total = (0, 0)
        for coefficient, z in zip(U, rows):
            total = add(total, mul((coefficient, 0), power(z, exponent)))
        assert total == (0, 0)

    raw = []
    for i, z in enumerate(rows):
        product = ONE
        for j, w in enumerate(rows):
            if i != j:
                product = mul(product, add(z, neg(w)))
        raw.append(div(power(z, 2), product))
    central, _ = primitive_vector(raw)
    multiplier = div(central[0], (U[0], 0))
    assert multiplier in ((1, 0), (-1, 0), (0, 1), (0, -1))
    assert central == [mul(multiplier, (u, 0)) for u in U]

    primitive, _ = primitive_vector(rows)
    radius_squared = sum(c*c for c in primitive[0])
    assert all(sum(c*c for c in z) == radius_squared for z in primitive)
    assert mul(primitive[0], primitive[1]) == mul(primitive[2], primitive[3])
    assert mul(primitive[0], primitive[1]) == mul(primitive[4], primitive[5])

    # Independent all-edge ordinary lcm verification of least radius.
    cotangents = []
    for z, w in combinations(rows, 2):
        product = mul(z, conj(w))
        cotangents.append(product[1]/(norms[0]-product[0]))
    scale = lcm(*(c.denominator for c in cotangents))
    qs = [int(c*scale) for c in cotangents]
    assert lcm(*(edge_norm(c, scale) for c in qs)) == radius_squared

    # Check the bad-prime error estimate as well as exact good-prime support.
    for a, b in ((2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
                 (5, 4), (7, 2), (6, 5), (8, 3), (8, 5), (9, 4)):
        prime = a*a+b*b
        exponents = [gaussian_valuation(z, (a, b)) for z in primitive]
        ordered = sorted(exponents)
        median_spread = sum(ordered[3:])-sum(ordered[:3])
        error_bound = 5*ordinary_valuation(2*scale, prime)
        for coefficient, exponent in zip(U, exponents):
            baseline = F(sum(abs(exponent-e) for e in exponents)-median_spread, 2)
            actual = ordinary_valuation(coefficient, prime)
            assert abs(actual-baseline) <= error_bound
    return q, radius_squared


def main():
    twice = ec_add(P, P)
    assert twice == (F(5476, 49), F(-135050, 343))
    assert (twice[0]-34).denominator != 1
    current = None
    small = []
    for multiple in range(1, 16):
        current = ec_add(current, P)
        assert current is not None
        q, n = configuration(current)
        if multiple % 2:
            assert 0 < current[0] < 3
            small.append((float(q), multiple))
        else:
            assert current[0] > 99
    best_q, best_multiple = min(small)
    print('PASS: 15 multiples, fixed primitive central coefficients, '
          'moments, Gaussian normalization, all-edge least radius, '
          'and central valuation errors at 11 split primes.')
    print(f'Smallest q among tested odd multiples: {best_q:.8g} at {best_multiple}P.')


if __name__ == '__main__':
    main()
