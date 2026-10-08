"""Audit the dual KZ sign directly on six-row ternary polynomials.

This finite certificate is independent of the published residue convention.
The general connection and coalition arguments are proved in
terminal_kernel_fusion_content.md.
"""

from fractions import Fraction
from itertools import permutations
from math import prod

from check_six_row_coherent_kernel_and_content import (
    CONVERSION, MATCHINGS, PARTITIONS,
)


def add(left, right, factor=1):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, 0) + factor * coefficient
        if not answer[monomial]:
            del answer[monomial]
    return answer


def permutation_sign(order):
    return (-1) ** sum(order[i] > order[j]
                      for i in range(3) for j in range(i + 1, 3))


def triangle_pair(partition):
    answer = {}
    for first in permutations(range(3)):
        for second in permutations(range(3)):
            assignment = [0] * 6
            for label, coordinate in zip(partition[0], first):
                assignment[label] = coordinate
            for label, coordinate in zip(partition[1], second):
                assignment[label] = coordinate
            answer[tuple(assignment)] = (
                permutation_sign(first) * permutation_sign(second))
    return answer


def source_basis():
    triangles = [triangle_pair(partition) for partition in PARTITIONS]
    basis = []
    for column in range(5):
        polynomial = {}
        for row in range(10):
            polynomial = add(polynomial, triangles[row],
                             CONVERSION[row][column])
        basis.append(polynomial)
    return basis


def swap_rows(polynomial, i, j):
    answer = {}
    for assignment, coefficient in polynomial.items():
        swapped = list(assignment)
        swapped[i], swapped[j] = swapped[j], swapped[i]
        answer[tuple(swapped)] = coefficient
    return answer


def audit(directions, basis):
    b = tuple(map(Fraction, directions))
    assert len(set(b)) == 6
    matching_values = [prod(b[j] - b[i] for i, j in matching)
                       for matching in MATCHINGS]
    polynomial = {}
    for value, generator in zip(matching_values, basis):
        polynomial = add(polynomial, generator, value)
    assert polynomial
    for i in range(6):
        derivative = {}
        for value, matching, generator in zip(matching_values, MATCHINGS, basis):
            u, v = next((u, v) for u, v in matching if i in (u, v))
            derivative = add(derivative, generator,
                             value * (-1 if i == u else 1) / (b[v] - b[u]))
        connection_term = {}
        for j in range(6):
            if j == i:
                continue
            casimir_pair = add(swap_rows(polynomial, i, j), polynomial,
                               -Fraction(1, 3))
            connection_term = add(connection_term, casimir_pair,
                                  1 / (4 * (b[i] - b[j])))
        scalar = sum((1 / (3 * (b[i] - b[j])) for j in range(6) if j != i),
                     Fraction(0))
        expected = add({}, polynomial, scalar)
        assert add(derivative, connection_term, -1) == expected
        # The opposite sign fails as a polynomial identity at these fixtures.
        assert add(derivative, connection_term, 1) != expected


def main():
    basis = source_basis()
    fixtures = [range(6), (5, 1, -2, 7, 0, 11),
                tuple(Fraction(i * i + i + 2, i + 2) for i in range(6))]
    for directions in fixtures:
        audit(directions, basis)
    print("PASS: dual KZ minus sign and pair-subtracted scalar in all 18 directions.")


if __name__ == "__main__":
    main()
