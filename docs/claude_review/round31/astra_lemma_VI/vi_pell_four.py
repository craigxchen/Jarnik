"""Verify the two four-point Pell scales by exact Gaussian arithmetic."""

from fractions import Fraction
from itertools import combinations
from math import prod


def multiply(z, w):
    x, y = z
    u, v = w
    return x*u - y*v, x*v + y*u


def norm(z):
    return z[0]**2 + z[1]**2


def conjugate(z):
    return z[0], -z[1]


def gaussian_gcd(z, w):
    while norm(w):
        numerator = multiply(z, conjugate(w))
        denominator = norm(w)
        q = tuple(round(Fraction(a, denominator)) for a in numerator)
        product = multiply(q, w)
        remainder = (z[0] - product[0], z[1] - product[1])
        assert norm(remainder) < denominator
        z, w = w, remainder
    return z


def advance(v, times):
    x, y = v
    for _ in range(times):
        x, y = 5*x + 2*y, 2*x + y
    return x, y


def four_point_test(rates, n, exponent):
    seeds = ((3, 1), (4, 1), (5, 1))
    blocks = [advance(v, rate*n) for v, rate in zip(seeds, rates)]
    points = []
    for flipped in ((), (0, 1), (0, 2), (1, 2)):
        point = (0, -1) if flipped else (1, 0)
        for j, z in enumerate(blocks):
            point = multiply(point, z if j in flipped else conjugate(z))
        points.append(point)
    raw_norm = prod(norm(z) for z in blocks)
    assert all(norm(z) == raw_norm for z in points)
    assert len(set(points)) == 4
    gcd = points[0]
    for z in points[1:]:
        gcd = gaussian_gcd(gcd, z)
    primitive_norm = raw_norm // norm(gcd)
    assert raw_norm == primitive_norm * norm(gcd)
    diameter_squared = max((x-u)**2 + (y-v)**2 for (x, y), (u, v) in combinations(points, 2))
    delta_squared = Fraction(diameter_squared, raw_norm)
    # Test the observed scale at these indices by rational comparisons.
    scale = delta_squared**exponent.denominator * primitive_norm**(2*exponent.numerator)
    assert Fraction(1, 2**40) < scale < 2**40
    return norm(gcd), primitive_norm.bit_length()


if __name__ == "__main__":
    for rates, exponent in [((1, 1, 1), Fraction(1, 3)), ((6, 7, 7), Fraction(3, 10))]:
        data = [four_point_test(rates, n, exponent) for n in (1, 2, 3, 4, 8, 16)]
        print("rates", rates, "Gaussian gcd norms and primitive norm bit lengths:", data)
