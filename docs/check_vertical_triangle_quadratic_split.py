#!/usr/bin/env python3
"""Exact bounded search on the vertical-line triangle surface.

This is finite evidence only. It uses Fraction arithmetic and tests whether
one real rational quadratic can split all seven source roots for each point
in the specified bounded family.
"""

from fractions import Fraction as F
from math import isqrt


def rational_values(numerator_bound, denominator_bound):
    return sorted(
        {
            F(a, b)
            for a in range(-numerator_bound, numerator_bound + 1)
            for b in range(1, denominator_bound + 1)
            if a != 0
        }
    )


def is_rational_square(x):
    return (
        x >= 0
        and isqrt(x.numerator) ** 2 == x.numerator
        and isqrt(x.denominator) ** 2 == x.denominator
    )


def is_gaussian_rational_square(x, y):
    """Whether x+i*y is a square in Q(i), using exact rational tests."""
    norm = x * x + y * y
    if not is_rational_square(norm):
        return False
    root_norm = F(isqrt(norm.numerator), isqrt(norm.denominator))
    return is_rational_square((root_norm + x) / 2) and is_rational_square(
        (root_norm - x) / 2
    )


def seven_roots(p, q, r, x):
    return [
        (F(0), F(1)),
        (x, p - 1),
        (x, q - 1),
        (x, r - 1),
        (x * (1 + 1 / p), p),
        (x * (1 + 1 / q), q),
        (x * (1 + 1 / r), r),
    ]


def main():
    triangle_values = rational_values(8, 4)
    triangles = []
    for p in triangle_values:
        for q in triangle_values:
            r = 1 - p - q
            if r == 0 or len({p, q, r}) < 3:
                continue
            if any(z in (0, 1, -1) for z in (p, q, r)):
                continue
            product = p * q * r
            if product <= 0 or not is_rational_square(product):
                continue
            x = F(isqrt(product.numerator), isqrt(product.denominator))
            triangles.append((p, q, r, x))

    # From i = d + K(ty + iy)^2, the common shift and leading
    # coefficient are d=(1/t-t)/2 and K=1/(2*t*y^2).
    t_values = rational_values(5, 4)
    y_values = rational_values(3, 3)
    candidate_count = 0
    for p, q, r, x in triangles:
        roots = seven_roots(p, q, r, x)
        for t in t_values:
            d = (1 / t - t) / 2
            for y in y_values:
                leading = 1 / (2 * t * y * y)
                candidate_count += 1
                if all(
                    is_gaussian_rational_square((real - d) / leading, imag / leading)
                    for real, imag in roots
                ):
                    print("FOUND", (p, q, r, x, d, leading))
                    return

    print(
        "No all-seven split in this finite search:",
        f"{len(triangles)} labeled surface points, {candidate_count} quadratics.",
    )


if __name__ == "__main__":
    main()
