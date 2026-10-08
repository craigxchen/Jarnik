#!/usr/bin/env python3
"""Compare relative-phase radius with a separately cleared full tuple."""

from itertools import product
from math import gcd, lcm

from check_positive_pentagon_primitive_content import (
    gaussian_div, gaussian_gcd, group_gcd, mul, norm, pentagon, product_ratio,
)


def conjugate(z):
    return z[0], -z[1]


def power(z, n):
    result = (1, 0)
    for _ in range(n):
        result = mul(result, z)
    return result


def reduced_phase(z):
    return gaussian_div(z, gaussian_gcd(z, conjugate(z)))


def gaussian_lcm(values):
    result = (1, 0)
    for z in values:
        result = gaussian_div(mul(result, z), gaussian_gcd(result, z))
    return result


def primitive_realization(rows):
    # Clear the ABSOLUTE phases, then remove the entire output gcd.
    scale = conjugate(gaussian_lcm([reduced_phase(z) for z in rows]))
    points = [gaussian_div(mul(scale, z), conjugate(z)) for z in rows]
    content = group_gcd(points)
    points = [gaussian_div(z, content) for z in points]
    assert norm(group_gcd(points)) == 1
    N = norm(points[0])
    assert all(norm(z) == N for z in points)
    return N, points


def check(rows):
    relatives = [mul(z, conjugate(rows[0])) for z in rows]
    relative_N = norm(gaussian_lcm([reduced_phase(z) for z in relatives]))
    N, points = primitive_realization(rows)
    assert relative_N == N
    pair_N = 1
    for i in range(5):
        for j in range(i):
            pair_N = lcm(pair_N, norm(reduced_phase(mul(rows[i], conjugate(rows[j])))))
    assert pair_N == N
    u0 = product_ratio(points, [(0, 3), (1, 2)], [(0, 2), (1, 3)])
    u1 = product_ratio(points, [(1, 4), (2, 3)], [(1, 3), (2, 4)])
    _, _, F = pentagon(u0.numerator, u0.denominator, u1.numerator, u1.denominator)
    assert F <= 2 * N
    return N, u0, u1, F


def main():
    # Include a half-angle pole, antipodal phases, unit phases, and odd-odd
    # rows explicitly; the power-of-(2+i) orbit fixture need not hit these.
    edge_rows = [(1, 0), (1, 1), (0, 1), (-1, 1), (-2, 1)]
    edge_base = check(edge_rows)
    assert check([(y, x) for x, y in edge_rows])[1:] == edge_base[1:]
    rows = [power((2, 1), j) for j in range(1, 6)]
    naive_N = norm(gaussian_lcm([reduced_phase(z) for z in rows]))
    base = check(rows)
    assert naive_N == 3125 and base[0] == 625
    explicit = [mul(power((2, 1), j), power((2, -1), 4 - j)) for j in range(5)]
    assert all(norm(z) == 625 for z in explicit)
    assert norm(group_gcd(explicit)) == 1
    checked = 0
    for a, b, c, d in product(range(-2, 3), repeat=4):
        if a * d == b * c or gcd(a, b, c, d) != 1:
            continue
        transformed = []
        for x, y in rows:
            X, Y = a * x + b * y, c * x + d * y
            content = gcd(X, Y)
            transformed.append((X // content, Y // content))
        result = check(transformed)
        assert result[1:] == base[1:]
        checked += 1
    print(f"PASS: {checked} projective matrices; relative lcm, full-tuple gcd, and all-edge lcm agree")
    print(f"Anchor counterexample: naive N={naive_N}, actual N={base[0]}, invariant F={base[3]}")


if __name__ == "__main__":
    main()
