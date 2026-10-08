#!/usr/bin/env python3
"""Exact finite audits of the two Pell endpoint-swapping parameters.

The infinite-family proof is in endpoint_swap_arithmetic_sensitivity.md.
Floating widths below are diagnostics; every asserted bound is rational.
"""

from itertools import combinations
from math import atan, exp, gcd, lcm, log

from check_fixed_pentagon_projective_radius import (
    conjugate, gaussian_div, gaussian_gcd, mul, norm, primitive_realization,
)
from check_integer_cotangent_normalization import check_clique


def all_edge_radius(rows):
    """Independent integer formula, with the odd-odd ramified correction."""
    radius = 1
    for (r, s), (u, v) in combinations(rows, 2):
        dot, cross = r * u + s * v, s * u - r * v
        content = gcd(dot, cross)
        a, b = dot // content, cross // content
        epsilon = 2 if a % 2 and b % 2 else 1
        radius = lcm(radius, (a * a + b * b) // epsilon)
    return radius


def full_radius(rows):
    radius, points = primitive_realization(rows)
    assert radius == all_edge_radius(rows)
    assert len(set(points)) == len(rows)
    return radius


def transform(rows, x, tau):
    images = [(x * r + tau * s, r - x * s) for r, s in rows]
    assert all(r > 0 and 0 <= x * s <= r for r, s in images)
    assert images[0] == (x, 1)
    assert images[1] == (x * x + tau, 0)
    for (r, s), (u, v) in zip(rows, images):
        assert (x * u + tau * v, u - x * v) == (
            (x * x + tau) * r, (x * x + tau) * s,
        )
    return images


def width(radius, x):
    return exp(log(radius) / 4 - log(x)) * (2 * x * atan(1 / x))


def check_alignment(rows, x, a, b, p, q):
    T, H = x * x + 1, max(p, q)
    images = [(q * x * r + (p * T - q * x * x) * s,
               q * r - q * x * s) for r, s in rows]
    assert all(r > 0 and 0 <= x * s <= r for r, s in images)
    assert images[0] == (q * x, q) and images[1] == (p * T, 0)
    F1, F2 = (q * x + 8 * p * a, q), (2 * q * x + 3 * p * b, 2 * q)
    assert images[2] == (b * F1[0], b * F1[1])
    assert images[3] == (a * F2[0], a * F2[1])
    assert q * x * (p * T - q * x * x) - q * q * x == q * x * (p - q) * T
    gaussian_div(mul(F1, F2), (x, 1))
    scale = conjugate(mul(F1, F2))
    common_points = [gaussian_div(mul(scale, z), conjugate(z)) for z in images]
    bound = norm(F1) * norm(F2)
    assert all(norm(z) == bound for z in common_points)
    target = full_radius(images)
    assert bound % target == 0 and target <= bound
    assert 8 * a < 13 * x and 8 * b < 5 * x
    assert norm(F1) < 197 * H * H * x * x
    assert 64 * norm(F2) < 1217 * H * H * x * x
    assert 239749 < 4 * 16**4
    assert 16 * target < (16 * H)**4 * x**4


def check_pell(index, U, V):
    assert U * U - 5 * V * V == 1 and U % 2 and V % 2 == 0
    a, b = 1 + 10 * V * V + 2 * U * V, 1 + 10 * V * V - 2 * U * V
    c = 13 + 130 * V * V + 38 * U * V
    x, u, w = 4 * U * V, 30 * U * V + 10 * V * V + 1, 20 * V * V + 16 * U * V + 2
    assert a - b == x and a * b == x * x + 1
    assert a * a - 3 * a * b + b * b == -1
    assert gcd(a, b) == 1 and a % 2 and b % 2
    rows = [(1, 0), (x, 1), (u, 8), (w, 3)]
    assert all(r > 0 and 0 <= x * s <= r for r, s in rows)
    xs = (3 * u, 24 * x, 8 * w)
    assert min(xs) == 24 * x
    source = full_radius(rows)
    assert source == check_clique(xs, 24) == 5 * a * b * c

    bad_rows = transform(rows, x, x * x)
    bad = full_radius(bad_rows)
    B0, B = bad_rows[0], bad_rows[2]
    assert b == u - 8 * x == 2 * U * (U - V) - 1
    assert b % V == 1 and b % U == U - 1 and gcd(b, 2 * x) == 1
    assert B == (x * (u + 8 * x), b)
    assert (B[0] - B0[0] * b, B[1] - B0[1] * b) == (16 * x * x, 0)
    assert gcd(*B) == 1 and B[0] % 2 == 0 and B[1] % 2
    for z in (B0, B):
        assert norm(gaussian_gcd(z, conjugate(z))) == 1
    assert norm(gaussian_gcd(B0, B)) == 1
    mandatory = norm(B0) * norm(B)
    assert bad % mandatory == 0
    assert mandatory >= x**4 * (u + 8 * x)**2
    assert u + 8 * x > 10 * V * V

    good_rows = transform(rows, x, x * x + 2)
    good = full_radius(good_rows)
    B1, B2 = (17 * a - b, 1), (a + 2 * b, 1)
    C2 = gaussian_div(B2, (1, 1))
    assert good_rows[2] == (b * B1[0], b * B1[1])
    assert good_rows[3] == (2 * a * B2[0], 2 * a * B2[1])
    assert good_rows[1] == (2 * a * b, 0)
    # B2/bar(B2) = i*C2/bar(C2), including the phase unit exactly.
    assert mul(B2, conjugate(C2)) == mul((0, 1), mul(C2, conjugate(B2)))
    for z in (B1, C2):
        assert norm(gaussian_gcd(z, conjugate(z))) == 1
    assert norm(gaussian_gcd(B1, C2)) == 1
    gaussian_div(mul(B1, C2), B0)
    d, f = 288 * a - 31 * b, (7 * a + 3 * b) // 2
    assert norm(B1) == a * d and norm(C2) == b * f
    assert good == a * b * d * f
    common_difference = 16 * a - 3 * b
    assert gcd(common_difference, b) == 1
    assert common_difference + 2 * f == 23 * a

    h, s = 1 + 10 * V * V, 2 * U * V
    assert h * h - 5 * s * s == 1 and s >= 72
    assert (x, a, b, d, f) == (2 * s, h + s, h - s, 257 * h + 319 * s, 5 * h + 2 * s)
    assert 16 * h * h < 81 * s * s
    assert 65 * 3589 * 53 < 256 * 15**4
    assert 16 * good < 15**4 * x**4  # Delta < 2/x gives C_good < 15.
    if index >= 11:
        assert 256 * source < x**4  # Source C < 1/2.
    if index == 1:
        assert source == 358853785
        assert mandatory == 2462387160981745
        assert bad == 22580168789266778196866305
    for p, q in ((1, 1), (1, 2), (2, 3), (2, 1), (3, 2), (3, 1), (5, 2), (7, 13)):
        check_alignment(rows, x, a, b, p, q)
    print(f"n={index}: C_source={width(source, x):.10g}, C_bad={width(bad, x):.10g}, C_good={width(good, x):.10g}")


def main():
    U, V, checked = 1, 0, 0
    for index in range(1, 32):
        U, V = 9 * U + 20 * V, 4 * U + 9 * V
        if index % 10 == 1:
            check_pell(index, U, V)
            checked += 1
    assert checked == 4
    print("PASS: four exact Pell cases; Gaussian/full-gcd and all-edge radii agree.")
    print("PASS: bad mandatory denominator, repaired exact radius, phase unit, and rational width bounds.")
    print("PASS: 32 rational-alignment cases, including negative matrix entries and the conformal control.")


if __name__ == "__main__":
    main()
