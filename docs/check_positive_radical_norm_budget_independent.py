#!/usr/bin/env python3
"""Independent exact audit of the positive-radical norm construction."""

from fractions import Fraction
from functools import reduce
from itertools import combinations, product
from math import comb, gcd, isqrt, lcm, prod


def rational_gcd(values):
    denominator = lcm(*(value.denominator for value in values))
    numerator = reduce(gcd, (abs(int(value * denominator)) for value in values))
    return Fraction(numerator, denominator)


def primitive_minors(rows):
    minors = {
        (i, j): rows[i][0] * rows[j][1] - rows[i][1] * rows[j][0]
        for i, j in combinations(range(len(rows)), 2)
    }
    denominator = lcm(*(value.denominator for value in minors.values()))
    integral = {pair: int(value * denominator) for pair, value in minors.items()}
    content = reduce(gcd, (abs(value) for value in integral.values()))
    return {pair: value // content for pair, value in integral.items()}


def bracket(minors, i, j):
    return minors[i, j] if i < j else -minors[j, i]


def check_exact_positive_normalization(point_count):
    # Pythagorean marked vectors make every radical rational, so the complete
    # sign and tau normalization can be checked with exact Fractions.
    vectors = [(u * u - 1, 2 * u) for u in range(2, point_count + 2)]
    determinant = lambda v, w: v[0] * w[1] - v[1] * w[0]
    values = [x * x + y * y for x, y in vectors]
    stars = [
        prod(determinant(vector, other) for j, other in enumerate(vectors) if i != j)
        for i, vector in enumerate(vectors)
    ]
    exponent = point_count - 2
    signed_rhos = [
        Fraction(isqrt(value) ** exponent, star) for value, star in zip(values, stars)
    ]
    rho_squares = [rho * rho for rho in signed_rhos]
    tau = rational_gcd(rho_squares)
    assert isqrt(tau.numerator) ** 2 == tau.numerator
    assert isqrt(tau.denominator) ** 2 == tau.denominator
    sqrt_tau = Fraction(isqrt(tau.numerator), isqrt(tau.denominator))
    radicands = [value / tau for value in rho_squares]
    assert all(value.denominator == 1 for value in radicands)
    radicands = [value.numerator for value in radicands]
    assert reduce(gcd, radicands) == 1
    assert all(isqrt(value) ** 2 == value for value in radicands)

    signs = [1 if star > 0 else -1 for star in stars]
    radical = sum(sign * isqrt(value) for sign, value in zip(signs, radicands))
    assert radical == sum(signed_rhos) / sqrt_tau > 0

    slopes = [Fraction(y, x) for x, y in vectors]
    divided_difference = sum(
        Fraction(isqrt(values[i]), vectors[i][0]) ** exponent
        / prod(slope - other for j, other in enumerate(slopes) if i != j)
        for i, slope in enumerate(slopes)
    )
    x_product = prod(x for x, _ in vectors)
    assert sum(signed_rhos) == divided_difference / x_product > 0

    if point_count == 7:
        # These are the actual two Gale columns.  Their zero sum is retained
        # simultaneously with the signs used above.
        gale_rows = [
            (
                Fraction(value**2 * vector[0], star),
                Fraction(value**2 * vector[1], star),
            )
            for value, star, vector in zip(values, stars, vectors)
        ]
        assert all(sum(row[column] for row in gale_rows) == 0 for column in (0, 1))
        minors = primitive_minors(gale_rows)
        assert reduce(gcd, (abs(value) for value in minors.values())) == 1
        square_stars = [
            -prod(bracket(minors, i, j) for j in range(7) if i != j)
            for i in range(7)
        ]
        assert all(value > 0 and isqrt(value) ** 2 == value for value in square_stars)
        roots = [isqrt(value) for value in square_stars]
        content = reduce(gcd, roots)
        assert [root // content for root in roots] == radicands


def extension_multiply(left, right, nonsquare, prime):
    return (
        (left[0] * right[0] + nonsquare * left[1] * right[1]) % prime,
        (left[0] * right[1] + left[1] * right[0]) % prime,
    )


def check_clean_cut_units():
    prime = 1009
    square_roots = {value * value % prime: value for value in range(prime)}
    root = square_roots[-1 % prime]
    nonsquare = next(value for value in range(2, prime) if value not in square_roots)

    def square_root(value):
        if value in square_roots:
            return square_roots[value], 0
        adjusted = value * pow(nonsquare, -1, prime) % prime
        return 0, square_roots[adjusted]

    for minority_size in (1, 2, 3):
        majority = list(range(7 - minority_size))
        scalings = []
        radicands = []
        for node in majority:
            derivative = prod(node - other for other in majority if other != node) % prime
            minority_factor = pow(node - root, minority_size, prime)
            full_derivative = derivative * minority_factor % prime
            q_value = (1 + node * node) % prime
            scaling = q_value**2 * pow(full_derivative, -1, prime) % prime
            scalings.append(scaling)
            radicands.append(
                pow(q_value, 5, prime) * pow(full_derivative, -2, prime) % prime
            )

        # A unit 2-by-2 minor proves that the limiting majority Gale rows
        # already span rank two modulo p, so saturation changes only units.
        determinant = scalings[0] * scalings[1] * (majority[1] - majority[0]) % prime
        assert determinant

        roots = [(0, 0)] * minority_size + [square_root(value) for value in radicands]
        signed_sum_product = (1, 0)
        for first_signs in product((-1, 1), repeat=6):
            signs = first_signs + (1,)
            value = (
                sum(sign * root_value[0] for sign, root_value in zip(signs, roots)) % prime,
                sum(sign * root_value[1] for sign, root_value in zip(signs, roots)) % prime,
            )
            assert value != (0, 0)
            signed_sum_product = extension_multiply(
                signed_sum_product,
                extension_multiply(value, value, nonsquare, prime),
                nonsquare,
                prime,
            )
        assert signed_sum_product[0] and signed_sum_product[1] == 0


def check_norm_budgets():
    seven_exponents = {degree: 38 * degree - 192 for degree in (1, 2, 4, 8, 16, 32, 64)}
    assert all(seven_exponents[degree] < 0 for degree in (1, 2, 4))
    assert seven_exponents[8] > 0

    for point_count in range(7, 102, 2):
        k = (point_count - 1) // 2
        a = sum(comb(2 * k, j) for j in range(k))
        r = 2 ** (point_count - 3)
        large = Fraction((point_count - 1) * r) - Fraction(point_count * a, 2)
        small = Fraction(point_count * a, 2)
        assert large > 0 and large + small == (point_count - 1) * r
        threshold = Fraction((point_count - 1) * r, large)
        possible_degrees = [2**power for power in range(point_count)]
        necessary = next(degree for degree in possible_degrees if degree >= threshold)
        assert 2 * large * necessary - 2 * (point_count - 1) * r >= 0
        if necessary > 1:
            assert 2 * large * (necessary // 2) - 2 * (point_count - 1) * r < 0
        # The maximum multiquadratic degree never contradicts the budget.
        assert 2 * large * 2 ** (point_count - 1) - 2 * (point_count - 1) * r > 0


def main():
    check_exact_positive_normalization(7)
    check_exact_positive_normalization(9)
    check_clean_cut_units()
    check_norm_budgets()
    print(
        "PASS: exact signs, tau normalization, Gale zero sum, clean-cut "
        "saturation, and odd-point norm budgets agree."
    )


if __name__ == "__main__":
    main()
