"""Exact checks for direct_primitive_axis_lift_obstruction.md."""

from fractions import Fraction
from itertools import product
from math import gcd

from check_gaussian_reflection_replacement import (
    gconj,
    gdivexact,
    ggcd,
    gmul,
    gnorm,
    gcd_many,
)


def gneg(x):
    return -x[0], -x[1]


def gprod(*values):
    out = (1, 0)
    for value in values:
        out = gmul(out, value)
    return out


def pell(index):
    u, v = 1, 0
    for _ in range(index):
        u, v = 9 * u + 20 * v, 4 * u + 9 * v
    assert u * u - 5 * v * v == 1
    return u, v


def family(index):
    u, v = pell(index)
    P = (-1, 2)
    A = (u + v, 2 * v)
    B = (2 * v, u - v)
    C = (2 * u + 8 * v, 3 * u + v)
    z = [
        gneg(gprod(P, A, B, C)),
        gneg(gprod(gconj(P), A, gconj(B), gconj(C))),
        gneg(gprod(gconj(P), gconj(A), B, gconj(C))),
        gneg(gprod(gconj(P), gconj(A), gconj(B), C)),
    ]
    return u, v, P, A, B, C, z


def real_after_multiplier(mu, z):
    return gmul(mu, z)[0]


def check_projection_bound(v, u, z):
    tested = 0
    for p, q in product(range(-30, 31), repeat=2):
        if p == q == 0:
            continue
        # The note writes mu=p-iq.
        mu = (p, -q)
        Q2 = p * p + q * q
        xs = [real_after_multiplier(mu, point) for point in z]
        width = max(xs) - min(xs)

        linear = (2 * p - q) * v - q * u
        assert width >= 16 * abs(linear)

        # Check w >= 16V/(5Q)-4Q/V without floating point.
        residual = Fraction(16 * v, 5) - Fraction(4 * Q2, v)
        if residual > 0:
            assert width * width * Q2 >= residual * residual
        if Q2 * 4 <= v * v:
            assert width * width * Q2 >= 4 * v * v
        tested += 1
    return tested


def check_clearing(anchor, gamma):
    common = ggcd(anchor, gamma)
    complement = gdivexact(anchor, common)
    d = gcd(abs(complement[0]), abs(complement[1]))
    h = gnorm(complement) // d
    quotient = gdivexact((h, 0), complement)
    multiplier = gdivexact(gmul((h, 0), gamma), anchor)
    assert gnorm(multiplier) == gnorm(quotient)
    return h, complement, multiplier


def check_index(index):
    u, v, P, A, B, C, z = family(index)
    assert index % 10 == 1 and v >= 4
    assert gnorm(gcd_many(z)) == 1
    assert all(gcd(abs(x), abs(y)) == 1 for x, y in z)

    a, b, c = gnorm(A), gnorm(B), gnorm(C)
    radius2 = 5 * a * b * c
    assert all(gnorm(point) == radius2 for point in z)

    expected_differences = [
        (0, 0),
        (32 * v, -16 * u - 16 * v),
        (-2 * u + 2 * v, 4 * v),
        (6 * u + 2 * v, -4 * u - 16 * v),
    ]
    differences = [
        (point[0] - z[0][0], point[1] - z[0][1]) for point in z
    ]
    assert differences == expected_differences

    tested = check_projection_bound(v, u, z)

    mu = gconj(B)
    relative_reals = [real_after_multiplier(mu, delta) for delta in differences]
    assert relative_reals == [0, -16, 0, -4]
    assert max(relative_reals) - min(relative_reals) == 16
    assert gnorm(mu) == b

    pair_data = [
        (0, gprod(P, B, C), A),
        (0, gprod(P, A, C), B),
        (0, gprod(P, A, B), C),
        (1, gprod(A, gconj(B)), gprod(gconj(P), gconj(C))),
        (1, gprod(A, gconj(C)), gprod(gconj(P), gconj(B))),
        (2, gprod(B, gconj(C)), gprod(gconj(P), gconj(A))),
    ]
    expected_imaginaries = [-8, 1, -1, -1, -3, 2]
    clearing = []
    for (anchor_index, gamma, expected_complement), imaginary in zip(
        pair_data, expected_imaginaries
    ):
        assert gamma[1] == imaginary
        assert gcd(abs(expected_complement[0]), abs(expected_complement[1])) == 1
        h, complement, multiplier = check_clearing(z[anchor_index], gamma)
        assert gnorm(complement) == gnorm(expected_complement)
        assert h == gnorm(expected_complement)
        assert gnorm(multiplier) == h
        clearing.append(h)

    gamma = gprod(P, A, C)
    assert gamma[1] == 1
    assert z[0] == gneg(gmul(B, gamma))
    h = b
    assert gdivexact(gmul((h, 0), gamma), z[0]) == gneg(gconj(B))
    lifted_radius2 = h * h * gnorm(gamma)
    assert lifted_radius2 == radius2 * b
    assert lifted_radius2 - (h * gamma[0]) ** 2 == b * b

    assert all(2 * v * v < height < 1200 * v * v for height in clearing)
    return tested


def check_generic_clearing_formula():
    cases = 0
    for x, y in product(range(-12, 13), repeat=2):
        if x == y == 0:
            continue
        d = gcd(abs(x), abs(y))
        h = (x * x + y * y) // d
        A = (x, y)
        gdivexact((h, 0), A)
        for smaller in range(1, h):
            numerator = gmul((smaller, 0), gconj(A))
            if numerator[0] % gnorm(A) == 0 and numerator[1] % gnorm(A) == 0:
                raise AssertionError((A, h, smaller))
        cases += 1
    return cases


def main():
    generic = check_generic_clearing_formula()
    projection = sum(check_index(index) for index in (1, 11, 21))
    print(f"PASS: {generic} exact least integer denominator-clearing cases.")
    print(f"PASS: {projection} exact integral multiplier projection bounds.")
    print("PASS: three primitive Pell tuples, coordinates, radii, and width-16 lifts.")
    print("PASS: all six pair-factor clearing scales and complementary contents.")
    print("PASS: PAC=X+i has cleared defect b^2 and squared radius R^2 b.")


if __name__ == "__main__":
    main()
