#!/usr/bin/env python3
"""Finite exact fixtures for real_algebraic_norm_descent_content.md.

These checks illustrate the prime-power resultant bound and collision
valuation vectors; the general assertions are proved in the note.
"""

from itertools import product
from math import gcd


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    return x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def conj(x):
    return x[0], -x[1]


def norm(x):
    return x[0] * x[0] + x[1] * x[1]


def nearest(num, den):
    q, r = divmod(num, den)
    return q + (2 * r > den)


def divmod_gaussian(x, y):
    den = norm(y)
    assert den
    real = x[0] * y[0] + x[1] * y[1]
    imag = x[1] * y[0] - x[0] * y[1]
    q = nearest(real, den), nearest(imag, den)
    return q, sub(x, mul(q, y))


def exact_quotient(x, y):
    q, r = divmod_gaussian(x, y)
    assert r == (0, 0), (x, y, r)
    return q


def gaussian_gcd(x, y):
    while y != (0, 0):
        _, r = divmod_gaussian(x, y)
        x, y = y, r
    return x


def product_gaussian(values):
    z = (1, 0)
    for value in values:
        z = mul(z, value)
    return z


def form(c, a, b):
    """A - i*c*B, whose conjugate is form(-c,A,B)."""
    return a, -c * b


def points(cs, a, b):
    vals = []
    for row in range(5):
        chosen = [
            form(-c if row and (mask & (1 << (row - 1))) else c, a, b)
            for mask, c in enumerate(cs, start=1)
        ]
        vals.append(product_gaussian(chosen))
    return vals


def exponents(cs):
    classes = sorted({c for c in cs} | {-c for c in cs})
    table = {c: [0] * 5 for c in classes}
    for row in range(5):
        for mask, c in enumerate(cs, start=1):
            selected = -c if row and (mask & (1 << (row - 1))) else c
            table[selected][row] += 1
    return table


def certificate(cs):
    table = exponents(cs)
    classes = sorted(table)
    mins = {c: min(table[c]) for c in classes}
    depth = max(table[c][row] - mins[c] for c in classes for row in range(5))
    resultant_product = 1
    for index, c in enumerate(classes):
        for other in classes[index + 1 :]:
            resultant_product *= abs(c - other) ** depth
    return table, mins, resultant_product


def check_fixture(name, cs, expected_vectors):
    table, mins, correction = certificate(cs)
    for c, vector in expected_vectors.items():
        actual = tuple(table[c][row] - table[c][0] for row in range(1, 5))
        assert actual == vector, (name, c, actual, vector)
    for c in table:
        if -c in table:
            sums = [table[c][row] + table[-c][row] for row in range(5)]
            assert len(set(sums)) == 1, (name, c, sums)

    checked = 0
    for a, b in product(range(-9, 10), range(-5, 6)):
        if gcd(a, b) != 1:
            continue
        vals = points(cs, a, b)
        G = (0, 0)
        for z in vals:
            G = gaussian_gcd(G, z)
        D = product_gaussian(
            form(c, a, b) for c in table for _ in range(mins[c])
        )
        residual = exact_quotient(G, D)
        exact_quotient((correction, 0), residual)

        for row in range(1, 5):
            Q = product_gaussian(
                form(c, a, b)
                for mask, c in enumerate(cs, start=1)
                if mask & (1 << (row - 1))
            )
            A = product_gaussian(
                form(c, a, b)
                for mask, c in enumerate(cs, start=1)
                if not mask & (1 << (row - 1))
            )
            assert sub(vals[row], vals[0]) == mul(A, sub(conj(Q), Q))
            assert norm(vals[row]) == norm(vals[0])
        checked += 1
    print(f"{name}: {checked} primitive parameter pairs, exact gcd and chords")


def main():
    distinct = list(range(1, 16))
    same = distinct.copy()
    same[1] = 1  # Labels {1} and {2} share the same oriented root.
    opposite = distinct.copy()
    opposite[1] = -1  # They share opposite orientations instead.

    check_fixture("distinct", distinct, {1: (-1, 0, 0, 0)})
    check_fixture("same-orientation collision", same, {1: (-1, -1, 0, 0)})
    check_fixture("opposite-orientation collision", opposite, {1: (-1, 1, 0, 0)})

    # N_{Q(sqrt(2),i)/Q(i)}(A+sqrt(2)*(t-i)) = A^2-2*(t-i)^2.
    # It has imaginary degree one although A is monic of degree eight.
    for t in range(-20, 21):
        A = (t * t + 1) * t**6
        square = mul((t, -1), (t, -1))
        descended = sub((A * A, 0), (2 * square[0], 2 * square[1]))
        assert descended[1] == 4 * t
    print("all finite exact fixtures passed")


if __name__ == "__main__":
    main()
