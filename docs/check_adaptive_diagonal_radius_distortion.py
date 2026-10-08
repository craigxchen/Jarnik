"""Exact checks for adaptive diagonal cotangent radius formulas."""

from itertools import combinations
from math import gcd


def lcm(a, b):
    return a // gcd(a, b) * b


def multiply(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conjugate(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def gaussian_gcd(z, w):
    while w != (0, 0):
        numerator = multiply(z, conjugate(w))
        denominator = norm(w)
        quotient = tuple((2 * x + denominator) // (2 * denominator)
                         for x in numerator)
        product = multiply(quotient, w)
        z, w = w, (z[0] - product[0], z[1] - product[1])
    return z


def divide_exact(z, w):
    numerator, denominator = multiply(z, conjugate(w)), norm(w)
    assert all(x % denominator == 0 for x in numerator)
    return tuple(x // denominator for x in numerator)


def independent_gaussian_radius(rows):
    # Construct the actual integral phases using Gaussian Euclidean lcms,
    # then remove their common Gaussian divisor. No edge-norm formula here.
    fractions = []
    common_denominator = (1, 0)
    for z in rows:
        common = gaussian_gcd(z, conjugate(z))
        numerator = divide_exact(z, common)
        denominator = divide_exact(conjugate(z), common)
        fractions.append((numerator, denominator))
        common_denominator = multiply(
            divide_exact(common_denominator,
                         gaussian_gcd(common_denominator, denominator)), denominator)
    points = [multiply(divide_exact(common_denominator, den), num)
              for num, den in fractions]
    content = points[0]
    for point in points[1:]:
        content = gaussian_gcd(content, point)
    primitive = [divide_exact(point, content) for point in points]
    assert len({norm(point) for point in primitive}) == 1
    return norm(primitive[0])


def edge_norm(num, den):
    g = gcd(abs(num), abs(den))
    num //= g
    den //= g
    eps = 2 if num % 2 and den % 2 else 1
    return (num * num + den * den) // eps


def direct_radius(xs, length, a, d):
    # First use the transformed finite cotangents and transformed residue.
    # The anchor is infinity, so its edge is represented separately.
    values = []
    for x in xs:
        values.append(edge_norm(a * x, d * length))
    for x, y in combinations(xs, 2):
        values.append(edge_norm(a * a * x * y + d * d * length * length,
                                a * d * length * (y - x)))
    out = 1
    for value in values:
        out = lcm(out, value)
    return out


def pair_formula(xs, length, a, d):
    values = []
    for x in xs:
        D, E = a * a * x, a * d * length
        g = gcd(abs(D), abs(E))
        D //= g
        E //= g
        eps = 2 if D % 2 and E % 2 else 1
        values.append((D * D + E * E) // eps)
    for x, y in combinations(xs, 2):
        D = a * a * x * y + d * d * length * length
        E = a * d * length * (y - x)
        g = gcd(abs(D), abs(E))
        D //= g
        E //= g
        eps = 2 if D % 2 and E % 2 else 1
        values.append((D * D + E * E) // eps)
    out = 1
    for value in values:
        out = lcm(out, value)
    return out


def radius_from_rows(rows, a, d):
    values = []
    for x, y in rows:
        values.append(edge_norm(a * x, d * y))
    for (x, y), (u, v) in combinations(rows, 2):
        values.append(edge_norm(a * a * x * u + d * d * y * v,
                                a * d * (x * v - y * u)))
    out = 1
    for value in values:
        out = lcm(out, value)
    return out


def main():
    cases = 0
    for length in range(1, 10):
        for xs in ((2, 5, 11), (4, 9, 16, 31), (6, 8, 16, 42)):
            for a, d in ((1, 1), (3, 1), (5, 2), (7, 3), (11, 4)):
                if gcd(a, d) != 1:
                    continue
                assert direct_radius(xs, length, a, d) == pair_formula(
                    xs, length, a, d
                )
                target = pair_formula(xs, length, a, d)
                transformed = [(a, 0)] + [(a * x, d * length) for x in xs]
                assert independent_gaussian_radius(transformed) == target
                for x, y in combinations(xs, 2):
                    dot = a * a * x * y + d * d * length * length
                    determinant = a * d * length * (y - x)
                    assert 2 * target * determinant ** 2 >= dot ** 2
                    assert target * (y - x) ** 2 >= 2 * x * y
                cases += 1

    matrix_cases = 0
    for rows in ([(1, 0), (2, 1), (3, 1), (-1, 2)],
                 [(1, 1), (0, 1), (-3, 1)],
                 [(4, 2), (12, -8), (3, 0)]):
        for matrix in (((2, 1), (1, 1)), ((1, 3), (-2, 1)),
                       ((0, 1), (1, 0)), ((3, 0), (0, 3)),
                       ((-1, 0), (0, 1))):
            transformed = [(matrix[0][0] * x + matrix[0][1] * y,
                            matrix[1][0] * x + matrix[1][1] * y)
                           for x, y in rows]
            predicted = 1
            for (x, y), (u, v) in combinations(transformed, 2):
                predicted = lcm(predicted, edge_norm(x * u + y * v, x * v - y * u))
            assert predicted == independent_gaussian_radius(transformed)
            matrix_cases += 1
    # The displayed Pell/double-inversion instance.
    assert pair_formula((3723, 3456, 7184), 24, 3, 1) == \
        2085984659911668125
    # Source reanchoring before a nonconformal map changes the map. For
    # H1=(2,1), the primitive transitions are (2,-1) and (7,-1).
    assert radius_from_rows([(1, 0), (2, 1), (3, 1)], 2, 1) == 629
    assert radius_from_rows([(1, 0), (2, -1), (7, -1)], 2, 1) == 3349
    print(f"PASS: {cases} diagonal lcm cases with independent Gaussian normalization "
          f"and exact pair bounds; {matrix_cases} general matrix cases; anchor correction.")


if __name__ == "__main__":
    main()
