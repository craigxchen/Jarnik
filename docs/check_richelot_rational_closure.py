"""Exact certificate for the Richelot rational-closure obstruction.

Only the Python standard library is used.  This tests a fixed rational
six-tuple, not an endpoint-profile construction.
"""

from fractions import Fraction
from itertools import combinations


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    value = Fraction(1)
    for c in range(len(a)):
        pivot = next((r for r in range(c, len(a)) if a[r][c]), None)
        if pivot is None:
            return 0
        if pivot != c:
            a[c], a[pivot] = a[pivot], a[c]
            value = -value
        pivot_value = a[c][c]
        value *= pivot_value
        for j in range(c, len(a)):
            a[c][j] /= pivot_value
        for r in range(c + 1, len(a)):
            factor = a[r][c]
            for j in range(c, len(a)):
                a[r][j] -= factor * a[c][j]
    assert value.denominator == 1
    return int(value)


def pairings(rows):
    if not rows:
        yield []
        return
    for j in range(1, len(rows)):
        for rest in pairings(rows[1:j] + rows[j + 1 :]):
            yield [(rows[0], rows[j])] + rest


def bracket(a, b):
    """Coefficients of a' b - a b', in increasing degree."""
    a0, a1, a2 = a
    b0, b1, b2 = b
    return [a1 * b0 - a0 * b1,
            2 * (a2 * b0 - a0 * b2),
            a2 * b1 - a1 * b2]


def discriminant(a):
    a0, a1, a2 = a
    return a1 * a1 - 4 * a0 * a2


def square_class(n):
    assert n != 0
    result = -1 if n < 0 else 1
    n = abs(n)
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent % 2:
            result *= p
        p += 1
    return result * n


def class_vector(n):
    """Coordinates in the squareclass span of -1,2,3,5."""
    bits = 1 if n < 0 else 0
    n = abs(n)
    for bit, prime in enumerate((2, 3, 5), 1):
        if n % prime == 0:
            bits |= 1 << bit
            n //= prime
    assert n == 1
    return bits


def binary_resultant(a, b):
    a0, a1, a2 = a
    b0, b1, b2 = b
    return determinant([[a2, a1, a0, 0], [0, a2, a1, a0],
                        [b2, b1, b0, 0], [0, b2, b1, b0]])


valid = 0
degenerate = 0
for pairs in pairings(list(range(6))):
    quadratics = [[a * b, -a - b, 1] for a, b in pairs]
    delta = determinant(quadratics)
    output = [bracket(quadratics[1], quadratics[2]),
              bracket(quadratics[2], quadratics[0]),
              bracket(quadratics[0], quadratics[1])]
    classes = [square_class(discriminant(q)) for q in output]
    if not delta:
        assert pairs == [(0, 5), (1, 4), (2, 3)]
        degenerate += 1
    else:
        assert all(binary_resultant(output[i], output[j])
                   for i, j in combinations(range(3), 2))
        vectors = [class_vector(c) for c in classes]
        # An ambient Galois sign choice fixing exactly two root pairs.
        signs = next((mask for mask in range(16)
                      if sum(bin(mask & v).count("1") % 2
                             for v in vectors) == 1), None)
        assert signs is not None
        valid += 1
    print(pairs, delta, classes)

assert (valid, degenerate) == (14, 1)
print("Verified: 14 nonsingular outputs each admit a one-pair Galois flip;")
print("the remaining pairing is degenerate.")
