"""Exact finite tests for scalar reduced chord reconstruction. No third-party dependencies."""

from fractions import Fraction
from itertools import combinations, product
from math import isqrt, prod


def separation(row, other, primes):
    return prod(p ** abs(a - b) for p, a, b in zip(primes, row, other))


def test_triangles():
    """Check the undivided triangle sine identity after two squarings."""
    primes = (5, 13, 17)
    rows = list(product(range(2), repeat=len(primes)))
    first = rows[0]
    tested = survivors = nonsquare = 0
    for second, third in combinations(rows[1:], 2):
        n = separation(first, second, primes)
        m = separation(first, third, primes)
        k = separation(second, third, primes)
        assert isqrt(n * m // k) ** 2 * k == n * m
        for a in range(1, isqrt(n) + 1):
            for b in range(1, isqrt(m) + 1):
                for c in range(1, isqrt(k) + 1):
                    tested += 1
                    x = n - a * a
                    y = m - b * b
                    difference = a * a * y + b * b * x - c * c * (n * m // k)
                    # Choosing either signed cosine admits both signs here.
                    if difference**2 != 4 * a * a * b * b * x * y:
                        continue
                    survivors += 1
                    nonsquare += int(isqrt(x) ** 2 != x or isqrt(y) ** 2 != y)
    assert survivors and nonsquare == 0
    print(f"Triangle magnitude triples tested: {tested}")
    print(f"Nonzero sine-law survivors: {survivors}")
    print(f"Survivors with a nonsquare real-component square: {nonsquare}")


def multiply(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def quotient_angle(z):
    a, b = z
    return Fraction(a * a - b * b, a * a + b * b), Fraction(2 * a * b, a * a + b * b)


def test_fan():
    primes = (5, 13, 17)
    bases = ((2, 1), (3, 2), (4, 1))
    profile = ((0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 1))
    valuations = ((0, 0, 0), (0, 1, 0), (1, 0, 0), (-1, 1, 1))
    w = []
    for row in valuations:
        number = (1, 0)
        for exponent, (a, b) in zip(row, bases):
            if exponent:
                number = multiply(number, (a, exponent * b))
        w.append(quotient_angle(number))
    assert len(set(w)) == len(w)
    for i, j in combinations(range(len(w)), 2):
        n = separation(profile[i], profile[j], primes)
        cos_angle = w[i][0] * w[j][0] + w[i][1] * w[j][1]
        square = n * (1 - cos_angle) / 2
        assert isqrt(square.numerator) ** 2 == square.numerator
        assert isqrt(square.denominator) ** 2 == square.denominator
        magnitude = Fraction(isqrt(square.numerator), isqrt(square.denominator))
        assert magnitude.denominator == (5 if (i, j) == (2, 3) else 1)
        print(f"Edge {i}{j}: separation norm {n:4d}, |scaled sine| {magnitude}")
    assert max(row[0] for row in valuations) - min(row[0] for row in valuations) == 2


if __name__ == "__main__":
    test_triangles()
    test_fan()
