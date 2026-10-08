#!/usr/bin/env python3
"""Exact radius/gcd audit on actual four-point source tuples.

The proof of asymptotic endpoint-cost inflation is in the adjacent note.
These fixtures do not claim an unbounded-point endpoint family.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[0], -a[1]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conjugate(a):
    return a[0], -a[1]


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def divide(a, b):
    numerator = mul(a, conjugate(b))
    return F(numerator[0], norm(b)), F(numerator[1], norm(b))


def exact_divide(a, b):
    q = divide(a, b)
    assert all(x.denominator == 1 for x in q)
    return tuple(x.numerator for x in q)


def gaussian_gcd(a, b):
    while b != (0, 0):
        q = divide(a, b)
        rounded = tuple((2 * x.numerator + x.denominator) // (2 * x.denominator) for x in q)
        a, b = b, sub(a, mul(rounded, b))
    return a


def gaussian_lcm(values):
    answer = (1, 0)
    for value in values:
        answer = exact_divide(mul(answer, value), gaussian_gcd(answer, value))
    return answer


def all_gcd(values):
    answer = (0, 0)
    for value in values:
        answer = gaussian_gcd(answer, value)
    return answer


def poly_mul(a, b):
    answer = [(0, 0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            answer[i + j] = add(answer[i + j], mul(x, y))
    return answer


def poly_add(a, b):
    return [add(a[i] if i < len(a) else (0, 0),
                b[i] if i < len(b) else (0, 0)) for i in range(max(len(a), len(b)))]


def evaluate(a, n):
    answer = (0, 0)
    for coefficient in reversed(a):
        answer = add(mul(answer, (n, 0)), coefficient)
    return answer


def determinant(a):
    answer = (1, 0)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j] != (0, 0)), None)
        if pivot is None:
            return 0, 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            answer = neg(answer)
        value = a[j][j]
        answer = mul(answer, value)
        for i in range(j + 1, len(a)):
            scale = divide(a[i][j], value)
            a[i] = [sub(x, mul(scale, y)) for x, y in zip(a[i], a[j])]
    return answer


def resultant(a, b):
    da, db = len(a) - 1, len(b) - 1
    rows = []
    for i in range(db):
        rows.append([(0, 0)] * i + list(reversed(a)) + [(0, 0)] * (db - 1 - i))
    for i in range(da):
        rows.append([(0, 0)] * i + list(reversed(b)) + [(0, 0)] * (da - 1 - i))
    return determinant(rows)


def primitive_factor(a, b, degree):
    numerator = a ** degree, b ** degree
    return exact_divide(numerator, (1, 1)) if a % 2 and b % 2 else numerator


def main():
    odd_power_cases = 0
    for a in range(-24, 25):
        for b in range(-24, 25):
            if gcd(a, b) != 1:
                continue
            source_factor = primitive_factor(a, b, 1)
            for degree in (1, 3, 5, 7, 9, 11):
                target_factor = primitive_factor(a, b, degree)
                divisor = (source_factor if degree % 4 == 1
                           else conjugate(source_factor))
                exact_divide(target_factor, divisor)
                assert norm(gaussian_gcd(target_factor,
                                         conjugate(target_factor))) == 1
                odd_power_cases += 1
    a_poly, b_poly = [(-1, 0), (1, 0), (1, 0)], [(1, 0), (2, 0)]
    q3 = poly_add(poly_add(poly_mul(a_poly, a_poly),
                          [mul((0, 1), z) for z in poly_mul(a_poly, b_poly)]),
                  [neg(z) for z in poly_mul(b_poly, b_poly)])
    polynomials = [[(-1, -1), (1, -2), (1, 0)],
                   [(-1, 0), (0, 1), (1, 0)],
                   [(0, 1), (2, 1), (1, 0)], q3]
    resultants = [resultant(a, b) for a, b in combinations(polynomials, 2)]
    assert resultants == [(6, 9), (6, -9), (225, 0), (-2, 0), (0, 9), (0, -9)]
    cost = 4264650
    norm_product = 1
    for value in resultants:
        norm_product *= norm(value)
    assert norm_product == cost ** 2
    e_poly = [(-1, 1), (1, 2), (1, 0)]
    improved_polynomials = [poly_mul(e_poly, e_poly),
                           [(0, -2), (1, 0)], [(1, -2), (1, 0)],
                           [(-1, -2), (1, -4), (1, 0)]]
    improved_resultants = [resultant(a, b) for a, b in combinations(improved_polynomials, 2)]
    assert improved_resultants == [(72, -54), (72, 54), (2025, 0), (1, 0), (3, 0), (3, 0)]
    improved_cost = 147622500
    norm_product = 1
    for value in improved_resultants:
        norm_product *= norm(value)
    assert norm_product == improved_cost ** 2

    for n in (20, 30, 40, 70, 100, 1000):
        signs = [(-1, -1), (1, -1), (-1, 1), (1, 1)]
        source = [exact_divide(mul((n, s), (n + 1, t)), (1, 1)) for s, t in signs]
        source_norm = ((n * n + 1) * ((n + 1) ** 2 + 1)) // 2
        assert all(norm(z) == source_norm for z in source)
        assert len(set(source)) == 4 and norm(all_gcd(source)) == 1
        aa, bb = n * n + n - 1, 2 * n + 1
        pairs = [(1, 0), (n, 1), (n + 1, 1), (aa, bb)]
        assert all(gcd(a, b) == 1 for a, b in pairs)
        for z, (a, b) in zip(source, pairs):
            assert mul(z, (a, -b)) == mul(source[0], (a, b))
        for degree in (1, 3, 5, 7):
            factors = [primitive_factor(a, b, degree) for a, b in pairs]
            assert all(norm(gaussian_gcd(h, conjugate(h))) == 1 for h in factors)
            anchor = gaussian_lcm([conjugate(h) for h in factors])
            output = [exact_divide(mul(anchor, (a ** degree, b ** degree)),
                                   (a ** degree, -b ** degree)) for a, b in pairs]
            assert norm(all_gcd(output)) == 1
            assert all(norm(z) == norm(anchor) for z in output)
            assert norm(anchor) >= source_norm
            if degree == 1:
                assert norm(anchor) == source_norm
            if degree != 3:
                continue
            c1, c2, c3 = factors[1:]
            cubic_norm = norm(anchor)
            numerator = norm(c1) * norm(c2) * norm(c3) * norm(all_gcd([c1, c2, c3]))
            denominator = 1
            for x, y in combinations([c1, c2, c3], 2):
                denominator *= norm(gaussian_gcd(x, y))
            assert numerator == cubic_norm * denominator
            values = [evaluate(p, n) for p in polynomials]
            assert c1 == mul((n, -1), values[1])
            assert c2 == exact_divide(mul((n + 1, -1), values[2]), (1, 1))
            assert c3 == exact_divide(mul(values[0], values[3]), (1, 1))
            for (x, y), res in zip(combinations(values, 2), resultants):
                exact_divide(res, gaussian_gcd(x, y))
            product_norm = 1
            for value in values:
                product_norm *= norm(value)
            assert product_norm <= 2 * cost ** 2 * cubic_norm
            assert cubic_norm <= product_norm
        if n % 20 == 0:
            factors = [(1, 0)]
            for a, b in pairs[1:]:
                real, imaginary = a ** 3 + 3 * a * b * b, 2 * b ** 3
                assert gcd(real, imaginary) == 2
                value = (real // 2, imaginary // 2)
                assert norm(gaussian_gcd(value, conjugate(value))) == 1
                assert mul(mul((a, b), (a, b)), (a, -2 * b)) == (real, imaginary)
                factors.append(value)
            anchor = gaussian_lcm([conjugate(h) for h in factors])
            output = [exact_divide(mul(anchor, value), conjugate(value)) for value in factors]
            assert norm(all_gcd(output)) == 1
            assert all(norm(z) == norm(anchor) for z in output)
            numerator = norm(all_gcd(factors[1:]))
            for value in factors[1:]:
                numerator *= norm(value)
            denominator = 1
            for x, y in combinations(factors[1:], 2):
                denominator *= norm(gaussian_gcd(x, y))
            assert numerator == norm(anchor) * denominator
            values = [evaluate(p, n) for p in improved_polynomials]
            product_norm = 1
            for value in values:
                product_norm *= norm(value)
            for (x, y), res in zip(combinations(values, 2), improved_resultants):
                exact_divide(res, gaussian_gcd(x, y))
            assert product_norm <= 4 * improved_cost ** 2 * norm(anchor)
            assert norm(anchor) <= product_norm
    assert 11 ** 4 < 2 * 10 ** 4
    print('PASS:', odd_power_cases, 'odd-power divisibilities, including negative slopes and antipodes.')
    print('PASS: six exact nonzero polynomial resultants and fixed gcd cost K=4264650.')
    print('PASS: six primitive four-point sources; exact half-angle ratios;')
    print('      primitive least-lcm outputs at d=1,3,5,7, no radius descent;')
    print('      complete cubic gcd formula and explicit two-sided radius bounds.')
    print('PASS: improved cubic 2u^3/(1+3u^2), six additional resultants,')
    print('      exact least-radius tuples and two-sided degree-eight bounds.')
    print('SCOPE: actual four-point endpoint counterfamily, not unbounded point count.')


if __name__ == '__main__':
    main()
