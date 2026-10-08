#!/usr/bin/env python3
"""Exact bookkeeping for terminal determinant degrees and height bounds."""
from fractions import Fraction
from math import comb, factorial

from check_coherent_nagata_balancing import (
    add_scaled, bounded_walk_count, clean, multiply, odd_determinant,
)


def evaluate_polynomial_pencil(poly_at_0, poly_at_1, parameter):
    result = dict(poly_at_0)
    add_scaled(result, poly_at_1, parameter)
    add_scaled(result, poly_at_0, -parameter)
    return clean(result)


def check_five_block_pencils():
    fixed = {i: Fraction(i * i + i + 2, i + 2) for i in range(5)}
    moving = 5
    fixed_polynomial = odd_determinant(tuple(range(5)), fixed)
    for anchor in range(5):
        block = (anchor,) + tuple(i for i in range(5) if i != anchor)
        # Permuting the determinant columns changes only this known sign.
        oriented = odd_determinant(block, fixed)
        mon = next(iter(fixed_polynomial))
        sign = Fraction(oriented[mon], fixed_polynomial[mon])
        assert sign in (1, -1)
        pencils = []
        for parameter in (Fraction(0), Fraction(1), Fraction(7, 3)):
            b = dict(fixed)
            b[moving] = parameter
            # Use the polynomial, denominator-cleared balancing coefficients.
            actual = {}
            for r, label in enumerate(block[1:], 1):
                rest = tuple(i for i in block[1:] if i != label)
                coefficient = Fraction((-1) ** (r + 1), sign)
                for other in rest:
                    coefficient *= b[other] - b[anchor]
                term = multiply(odd_determinant((anchor, moving, label), b),
                                odd_determinant(rest, b))
                add_scaled(actual, term, coefficient)
            expected = {}
            add_scaled(expected, fixed_polynomial, parameter - b[anchor])
            assert clean(actual) == clean(expected)
            pencils.append(clean(actual))
        assert evaluate_polynomial_pencil(pencils[0], pencils[1], Fraction(7, 3)) == pencils[2]


def main():
    check_five_block_pencils()
    table = []
    for k in range(1, 41):
        d, m = 2 * k, 6 * k
        n = m - 1
        rho = bounded_walk_count(m, d)
        delta = bounded_walk_count(n, d)
        beta = bounded_walk_count(n, d - 1)
        assert rho == delta + beta
        assert m * delta + m * beta - rho * d == 2 * d * rho
        assert rho == comb(m, d) - 2 * comb(m, d - 1) + comb(m, d - 2)
        assert delta == comb(n, d) - 2 * comb(n, d - 1) + comb(n, d - 2)
        assert delta == Fraction(comb(m, d) * (d - 1), 3 * (2 * d + 1))
        r = 2 * factorial(m) // (factorial(d) * factorial(d + 1) * factorial(d + 2))
        assert r >= rho + delta
        majority = sum(comb(m, t) * max(0, t - m // 2) for t in range(m + 1))
        assert majority == Fraction(m * comb(m, m // 2), 4)
        raw = delta * m * 2 ** (m - 3)
        content = delta * majority
        A0 = m * 2 ** (m - 2) - (m // 2) * comb(m, m // 2)
        assert raw - content == Fraction(delta * A0, 2)
        height = Fraction(delta * A0, 2 * (r - rho))
        ratio = height / 2 ** (d + 1)
        assert ratio > 1
        if k >= 2:
            lower = Fraction(d * (d + 1)**2 * (d + 2) * (d - 1), 64 * (2*d + 1))
            assert ratio > lower >= Fraction(25, 8)
        if k <= 5:
            table.append((k, delta, rho, r - rho, float(height), float(ratio)))
    print('PASS: moving-five-block pencil identity for all five anchors')
    print('PASS: image/source degree split, majority content, and height formulas k=1..40')
    print('k, delta, rho, kernel rank, first-minimum exponent, threshold ratio:')
    for row in table:
        print(row)


if __name__ == '__main__':
    main()
