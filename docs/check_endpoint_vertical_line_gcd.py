"""Exact checks for endpoint_vertical_line_gcd.md."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt


Gaussian = tuple[int, int]


def add(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] + w[0], z[1] + w[1]


def sub(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] - w[0], z[1] - w[1]


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def norm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def divmod_gaussian(z: Gaussian, w: Gaussian) -> tuple[Gaussian, Gaussian]:
    denominator = norm(w)
    q_real = round(Fraction(z[0] * w[0] + z[1] * w[1], denominator))
    q_imag = round(Fraction(z[1] * w[0] - z[0] * w[1], denominator))
    quotient = q_real, q_imag
    return quotient, sub(z, mul(quotient, w))


def exact_div(z: Gaussian, w: Gaussian) -> Gaussian:
    quotient, remainder = divmod_gaussian(z, w)
    assert remainder == (0, 0), (z, w, quotient, remainder)
    return quotient


def gcd_gaussian(z: Gaussian, w: Gaussian) -> Gaussian:
    while w != (0, 0):
        _, remainder = divmod_gaussian(z, w)
        z, w = w, remainder
    return z


def lcm_gaussian(z: Gaussian, w: Gaussian) -> Gaussian:
    if z == (0, 0) or w == (0, 0):
        return 0, 0
    return mul(exact_div(z, gcd_gaussian(z, w)), w)


def gcd_many(values: list[Gaussian]) -> Gaussian:
    result = (0, 0)
    for value in values:
        result = gcd_gaussian(result, value)
    return result


def lcm_many(values: list[Gaussian]) -> Gaussian:
    result = (1, 0)
    for value in values:
        result = lcm_gaussian(result, value)
    return result


def check_pair_formula(limit: int = 30) -> int:
    checks = 0
    for d in range(1, limit + 1):
        for t_i in range(-limit, limit + 1):
            for t_j in range(t_i + 1, limit + 1):
                delta = t_j - t_i
                w_i = -d, t_i
                w_j = -d, t_j
                actual = norm(gcd_gaussian(w_i, w_j))
                predicted = gcd(
                    gcd(delta * delta, abs(delta * d)),
                    gcd(abs(delta * t_i), d * d + t_i * t_i),
                )
                assert actual == predicted, (d, t_i, t_j, actual, predicted)
                assert actual <= abs(delta) * min(d, abs(delta))
                checks += 1
    return checks


def circle_points(squared_radius: int) -> list[Gaussian]:
    bound = isqrt(squared_radius)
    points = []
    for x in range(-bound, bound + 1):
        y_squared = squared_radius - x * x
        y = isqrt(y_squared)
        if y * y != y_squared:
            continue
        points.append((x, y))
        if y:
            points.append((x, -y))
    return points


def check_full_tuple_reconstruction(limit: int = 500) -> int:
    checks = 0
    for squared_radius in range(1, limit + 1):
        points = circle_points(squared_radius)
        if len(points) < 2 or norm(gcd_many(points)) != 1:
            continue
        for anchor in points:
            twice_anchor = mul((2, 0), anchor)
            cofactors: list[Gaussian] = []
            primitive_chords: list[Gaussian] = []
            for point in points:
                if point == anchor:
                    continue
                chord = sub(point, anchor)
                content = gcd(abs(chord[0]), abs(chord[1]))
                primitive = chord[0] // content, chord[1] // content
                cofactor = exact_div(twice_anchor, primitive)
                assert cofactor[0] == -content
                cofactors.append(cofactor)
                primitive_chords.append(primitive)

            cofactor_lcm = lcm_many(cofactors)
            chord_gcd = gcd_many(primitive_chords)
            quotient = exact_div(twice_anchor, cofactor_lcm)
            assert norm(quotient) == norm(chord_gcd)
            assert norm(chord_gcd) <= norm((2, 0))
            assert squared_radius <= norm(cofactor_lcm) <= 4 * squared_radius

            for left, right in combinations(cofactors, 2):
                if left[0] != right[0]:
                    continue
                d = -left[0]
                delta = right[1] - left[1]
                if delta == 0:
                    continue
                actual = norm(gcd_gaussian(left, right))
                predicted = gcd(
                    gcd(delta * delta, abs(delta * d)),
                    gcd(abs(delta * left[1]), d * d + left[1] * left[1]),
                )
                assert actual == predicted
            checks += 1
    return checks


if __name__ == "__main__":
    pair_checks = check_pair_formula()
    tuple_checks = check_full_tuple_reconstruction()
    print(f"PASS: {pair_checks} vertical-line gcd checks")
    print(f"PASS: {tuple_checks} primitive full-circle reconstruction checks")
