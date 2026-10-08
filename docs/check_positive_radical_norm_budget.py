"""Exact checks of the positive-radical norm-budget notes; standard library only."""

from fractions import Fraction
from itertools import product
from math import comb


def trim(poly):
    while len(poly) > 1 and not poly[-1]:
        poly.pop()
    return poly


def add(left, right):
    out = [Fraction(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] += value
    return trim(out)


def multiply(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def derivative(poly):
    return [i * value for i, value in enumerate(poly)][1:] or [Fraction(0)]


def check_derivatives():
    q = list(map(Fraction, (5, 2, 3)))
    determinant = q[0] * q[2] - q[1] ** 2 / 4
    for n in range(3, 32, 2):
        exponent = Fraction(n - 2, 2)
        numerator = [Fraction(1)]
        for order in range(n - 1):
            numerator = add(
                multiply(q, derivative(numerator)),
                [
                    (exponent - order) * value
                    for value in multiply(derivative(q), numerator)
                ],
            )
        double_factorial = 1
        for odd in range(1, n - 1, 2):
            double_factorial *= odd
        assert numerator == [double_factorial ** 2 * determinant ** ((n - 1) // 2)]


def check_fair_counts():
    for n in range(7, 102, 2):
        k = (n - 1) // 2
        a = sum(comb(n - 1, j) for j in range(k))
        b = sum(comb(n - 2, j) for j in range(k - 1))
        r = 2 ** (n - 3)
        assert a == (4 ** k - comb(2 * k, k)) // 2
        assert a - b == r
        large = Fraction((n - 2) * a, 2) - (n - 1) * b
        tiny = Fraction(n * a, 2)
        assert large + tiny == (n - 1) * r
        assert large > 0
        assert large * (2 ** (n - 1) - 1) - tiny > 0
        for minority in range(1, k + 1):
            assert n - 2 - 2 * (minority - 1) == n - 2 * minority > 0


def check_clean_units():
    prime = 1009
    square_roots = {x * x % prime: x for x in range(prime)}
    root_minus_one = square_roots[prime - 1]
    nonsquare = next(x for x in range(2, prime) if x not in square_roots)

    def multiply_extension(z, w):
        return (
            (z[0] * w[0] + nonsquare * z[1] * w[1]) % prime,
            (z[0] * w[1] + z[1] * w[0]) % prime,
        )

    def square_root(value):
        if value in square_roots:
            return (square_roots[value], 0)
        return (0, square_roots[value * pow(nonsquare, -1, prime) % prime])

    expected = {
        1: ([777, 15, 259, 806, 457, 795], 74),
        2: ([755, 224, 201, 506, 801], 227),
        3: ([28, 469, 829, 977], 951),
    }
    assert root_minus_one == 540
    for minority in (1, 2, 3):
        nodes = list(range(7 - minority))
        radicands = []
        for node in nodes:
            denominator = pow(node - root_minus_one, minority, prime)
            for other in nodes:
                if node != other:
                    denominator = denominator * (node - other) % prime
            radicands.append(
                pow(1 + node * node, 5, prime) * pow(denominator, -2, prime) % prime
            )
        assert radicands == expected[minority][0]
        roots = [(0, 0)] * minority + [square_root(value) for value in radicands]
        norm = (1, 0)
        for first_signs in product((-1, 1), repeat=6):
            signs = first_signs + (1,)
            value = (
                sum(sign * root[0] for sign, root in zip(signs, roots)) % prime,
                sum(sign * root[1] for sign, root in zip(signs, roots)) % prime,
            )
            assert value != (0, 0)
            norm = multiply_extension(norm, multiply_extension(value, value))
        assert norm == (expected[minority][1], 0)


def check_determinant_conjugates():
    cut_counts = {1: 7, 2: 21, 3: 35}
    expected_content = {1: Fraction(0), 2: Fraction(-35, 2), 3: Fraction(-63)}
    expected_totals = {1: Fraction(1120), 2: Fraction(896), 3: Fraction(0)}
    sign_counts = [1, 7, 21, 35]
    for size in (1, 2, 3):
        content = Fraction(0)
        for minority, count in cut_counts.items():
            minimum = min(
                a * (a + Fraction(7, 2) - minority - size)
                for a in range(max(0, size - 7 + minority), min(size, minority) + 1)
            )
            content += count * minimum
        assert content == expected_content[size]
        total = Fraction(0)
        for flips, count in enumerate(sign_counts):
            selected = min(size, flips)
            raw = -77 * (size - selected) + 19 * selected - 16 * selected * (selected - 1)
            total += count * (raw - content)
        assert total == expected_totals[size]

    def bit_count(value):
        return bin(value).count("1")

    exponents = {}
    for bits in range(64):
        flips = min(bit_count(bits), 7 - bit_count(bits))
        signs = [-1 if bits & (1 << i) else 1 for i in range(7)]
        value = -168 + 112 * flips - 16 * flips * flips
        assert value == -8 * sum(signs[i] * signs[j] for i in range(7) for j in range(i + 1, 7))
        exponents[bits] = value
    for character in range(64):
        weight = bit_count(character)
        weight += weight % 2  # Include label 7 to make the character even.
        coefficient = Fraction(
            sum(value * (-1 if bit_count(bits & character) % 2 else 1)
                for bits, value in exponents.items()),
            64,
        )
        assert coefficient == (-8 if weight == 2 else 0)


if __name__ == "__main__":
    check_derivatives()
    check_fair_counts()
    check_clean_units()
    check_determinant_conjugates()
    print("PASS: exact positive derivative identities for odd n=3,...,31.")
    print("PASS: fair incidence and complete conjugate budgets for odd n=7,...,101.")
    print("PASS: all 64 squared signed sums are units at each of three clean cut types modulo 1009.")
    print("PASS: determinant content shifts, complete sign-class budgets, and all 64 pair-character Fourier coefficients.")
