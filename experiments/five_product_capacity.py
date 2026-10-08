#!/usr/bin/env python3
"""Exact certificates for five-product reduced-determinant capacity.

The ten variables ``x_ij`` are nonnegative multiplicities of reduced pair
determinants.  An atom with signed occurrence word ``eps`` charges pair ``ij``
by ``abs(eps[i] - eps[j]) / 2``.  A one-use conductor budget imposes one such
linear constraint for every occurrence word.

No floating-point linear-programming package is used.  The displayed primal
and dual vectors are checked with ``fractions.Fraction``; equality of their
objectives proves optimality.
"""

from fractions import Fraction as Q
from itertools import combinations, product
from math import atan2, gcd, log


N = 5
EDGES = tuple(combinations(range(N), 2))


def gaussian_mul(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def gaussian_conj(z):
    return z[0], -z[1]


def det(matrix):
    """Exact determinant by fraction-free Laplace expansion (tiny matrices)."""
    n = len(matrix)
    assert all(len(row) == n for row in matrix)
    if n == 0:
        return 1
    if n == 1:
        return matrix[0][0]
    return sum((-1) ** j * matrix[0][j]
               * det(tuple(row[:j] + row[j + 1:] for row in matrix[1:]))
               for j in range(n))


def cost(eps, edge):
    i, j = edge
    return Q(abs(eps[i] - eps[j]), 2)


def all_ternary():
    return tuple(product((-1, 0, 1), repeat=N))


def rooted_binary():
    return tuple(product((0, 1), repeat=N))


def chronological_alternating():
    words = [(0,) * N]
    for support_bits in product((0, 1), repeat=N):
        support = [i for i, bit in enumerate(support_bits) if bit]
        if not support:
            continue
        for first in (-1, 1):
            eps = [0] * N
            for k, i in enumerate(support):
                eps[i] = first * (-1) ** k
            words.append(tuple(eps))
    assert len(words) == len(set(words)) == 63
    return tuple(words)


def check_primal(patterns, x, claimed):
    assert len(x) == len(EDGES)
    assert all(weight >= 0 for weight in x)
    assert sum(x) == claimed
    for eps in patterns:
        assert sum(weight * cost(eps, edge)
                   for weight, edge in zip(x, EDGES)) <= 1


def check_dual(patterns, y, claimed):
    pattern_set = set(patterns)
    assert all(eps in pattern_set and weight >= 0
               for eps, weight in y.items())
    assert sum(y.values()) == claimed
    for edge in EDGES:
        assert sum(weight * cost(eps, edge)
                   for eps, weight in y.items()) >= 1


def literal_integer_optimum(patterns, edge_cap):
    """Brute-force the small literal-multiplicity integer program."""
    best_total = -1
    best = None
    for multiplicities in product(range(edge_cap + 1), repeat=len(EDGES)):
        total = sum(multiplicities)
        if total <= best_total:
            continue
        if all(sum(Q(m) * cost(eps, edge)
                   for m, edge in zip(multiplicities, EDGES)) <= 1
               for eps in patterns):
            best_total = total
            best = multiplicities
    return best_total, best


def ternary_certificate():
    patterns = all_ternary()
    x = (Q(1, 6),) * 10
    y = {}
    for plus_tuple in combinations(range(N), 2):
        plus = set(plus_tuple)
        eps = tuple(1 if i in plus else -1 for i in range(N))
        y[eps] = Q(1, 6)
    check_primal(patterns, x, Q(5, 3))
    check_dual(patterns, y, Q(5, 3))
    return x, y


def binary_certificate():
    patterns = rooted_binary()
    x = (Q(1, 3),) * 10
    y = {}
    for ones_tuple in combinations(range(N), 2):
        ones = set(ones_tuple)
        eps = tuple(1 if i in ones else 0 for i in range(N))
        y[eps] = Q(1, 3)
    check_primal(patterns, x, Q(10, 3))
    check_dual(patterns, y, Q(10, 3))
    return x, y


def alternating_certificate():
    patterns = chronological_alternating()
    # Edge order: 01,02,03,04,12,13,14,23,24,34.
    x = (Q(0), Q(1, 2), Q(1, 2), Q(0), Q(1, 2),
         Q(1, 2), Q(0), Q(0), Q(0), Q(0))
    y = {
        (-1, 0, 0, 1, -1): Q(1, 2),
        (-1, 0, 1, -1, 1): Q(1, 2),
        (-1, 1, 0, 0, -1): Q(1, 2),
        (-1, 1, -1, 0, 1): Q(1, 2),
    }
    check_primal(patterns, x, Q(2))
    check_dual(patterns, y, Q(2))
    return x, y


def rank_four_saturated_counter_support():
    support = []
    for plus_tuple in combinations(range(N), 2):
        plus = set(plus_tuple)
        support.append(tuple(1 if i in plus else -1 for i in range(N)))
    for r in range(1, N):
        support.append(tuple(1 if i == r else 0 for i in range(N)))

    differences = [tuple(eps[r] - eps[0] for r in range(1, N))
                   for eps in support]
    for r in range(4):
        target = tuple(1 if j == r else 0 for j in range(4))
        assert target in differences
    assert all(not (len(set(eps)) == 1 and eps[0] != 0)
               for eps in support)
    return tuple(support)


def squarefree_six_point_certificate():
    """Exact 13/6 certificate for the squarefree six-point calibration."""
    primes = (5, 13, 17, 29, 41, 61, 73, 89)
    gaussian_primes = (
        (2, 1), (3, 2), (4, 1), (5, 2),
        (5, 4), (6, 5), (8, 3), (8, 5),
    )
    transition_rows = (
        (-1, -1, -1, +1, 0, 0, +1, +1),
        (0, 0, +1, -1, 0, +1, 0, -1),
        (+1, +1, -1, 0, +1, -1, 0, 0),
        (-1, -1, +1, +1, 0, 0, -1, 0),
        (0, +1, 0, -1, -1, +1, 0, +1),
    )
    patterns = tuple(zip(*transition_rows))

    # Edge order: 01,02,03,04,12,13,14,23,24,34.
    x = (Q(0), Q(1, 2), Q(1, 2), Q(1, 3), Q(0),
         Q(1, 6), Q(1, 6), Q(0), Q(1, 2), Q(0))
    y_values = (Q(1, 3), Q(1, 12), Q(7, 12), Q(5, 12),
                Q(0), Q(0), Q(1, 12), Q(2, 3))
    y = dict(zip(patterns, y_values))
    check_primal(patterns, x, Q(13, 6))
    check_dual(patterns, y, Q(13, 6))

    conductor = 1
    for p in primes:
        conductor *= p
    assert conductor == 520_699_108_865
    points = (
        (707_201, 143_408),
        (706_841, 145_172),
        (706_663, 146_036),
        (706_288, 147_839),
        (706_151, 148_492),
        (705_944, 149_473),
    )
    assert all(a * a + b * b == conductor and gcd(a, b) == 1
               for a, b in points)

    raw_products = []
    for row in transition_rows:
        value = (1, 0)
        for state, prime in zip(row, gaussian_primes):
            if state == 1:
                value = gaussian_mul(value, prime)
            elif state == -1:
                value = gaussian_mul(value, gaussian_conj(prime))
        raw_products.append(value)
    assert tuple(raw_products) == (
        (14_429, 18),
        (1_636, 1),
        (1_174, 1_177),
        (1_082, -1_081),
        (6_473, 6_482),
    )
    half_gap_products = (
        raw_products[0],
        raw_products[1],
        gaussian_mul(raw_products[2], (1, -1)),
        gaussian_mul(raw_products[3], (1, 1)),
        gaussian_mul(raw_products[4], (1, -1)),
    )
    assert half_gap_products == (
        (14_429, 18), (1_636, 1), (2_351, 3),
        (2_163, 1), (12_955, 9),
    )

    reduced_adjacent = []
    for (a, b), (c, d) in zip(points, points[1:]):
        dot = a * c + b * d
        wedge = a * d - b * c
        content = gcd(dot, wedge)
        reduced_adjacent.append((dot // content, wedge // content))
    expected_squares = []
    for index, value in enumerate(half_gap_products):
        square = gaussian_mul(value, value)
        divisor = 1 if index < 2 else 2
        expected_squares.append((square[0] // divisor, square[1] // divisor))
    assert reduced_adjacent == expected_squares

    difference_minor = tuple(
        tuple(transition_rows[r][c] - transition_rows[0][c]
              for c in range(4))
        for r in range(1, 5)
    )
    assert det(difference_minor) == -6

    arguments = [atan2(b, a) for a, b in half_gap_products]
    diameter = max(arguments) - min(arguments)
    coefficient_four_slack = log(conductor) - 4 * log(1 / diameter)
    assert abs(diameter - 0.000813731233179208) < 1e-15
    assert abs(coefficient_four_slack + 1.4770835217147) < 1e-12
    return x, y


def main():
    certificates = (
        ("ternary", Q(5, 3), ternary_certificate()),
        ("binary", Q(10, 3), binary_certificate()),
        ("chronological", Q(2), alternating_certificate()),
    )
    for name, optimum, (x, y) in certificates:
        print(name, "fractional optimum", optimum,
              "primal", x, "dual support", len(y))

    x, y = squarefree_six_point_certificate()
    print("squarefree calibration fractional optimum", Q(13, 6),
          "primal", x, "dual support", len(y))

    integer_models = (
        ("ternary", all_ternary(), 1),
        ("binary", rooted_binary(), 2),
        ("chronological", chronological_alternating(), 1),
    )
    for name, model, cap in integer_models:
        optimum, witness = literal_integer_optimum(model, cap)
        print(name, "literal integer optimum", optimum, "witness", witness)

    support = rank_four_saturated_counter_support()
    print("rank-four saturated counter-support words", len(support))


if __name__ == "__main__":
    main()
