#!/usr/bin/env python3
"""Exact checks for the eight-point quadratic conic-boundary barrier."""

from fractions import Fraction
from itertools import combinations
from math import gcd


TRIPLES = list(combinations(range(8), 3))


def det3(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def evaluations(monomials, cut, parameters):
    rows = [
        (1, parameters[i], 0) if i in cut else (1, 0, parameters[i])
        for i in range(8)
    ]
    needed = {triple for pair in monomials for triple in pair}
    pluckers = {
        triple: det3([rows[i] for i in triple]) for triple in needed
    }
    return [pluckers[left] * pluckers[right] for left, right in monomials]


def primitive_row(row):
    divisor = 0
    for value in row:
        divisor = gcd(divisor, abs(value))
    return tuple(value // divisor for value in row)


def determinant(matrix):
    matrix = [list(map(Fraction, row)) for row in matrix]
    answer = Fraction(1)
    for column in range(len(matrix)):
        pivot = next(
            row
            for row in range(column, len(matrix))
            if matrix[row][column]
        )
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            answer = -answer
        value = matrix[column][column]
        answer *= value
        for row in range(column + 1, len(matrix)):
            multiplier = matrix[row][column] / value
            for j in range(column, len(matrix)):
                matrix[row][j] -= multiplier * matrix[column][j]
    return answer


def add_polynomials(left, right, scale=1):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, 0) + scale * coefficient
        if answer[monomial] == 0:
            del answer[monomial]
    return answer


def multiply_polynomials(left, right):
    answer = {}
    for powers_left, coefficient_left in left.items():
        for powers_right, coefficient_right in right.items():
            powers = tuple(a + b for a, b in zip(powers_left, powers_right))
            answer[powers] = answer.get(powers, 0) + coefficient_left * coefficient_right
    return {powers: coefficient for powers, coefficient in answer.items() if coefficient}


def bracket(i, j, variables):
    positive = [0] * variables
    negative = [0] * variables
    positive[j] = 1
    negative[i] = 1
    return {tuple(positive): 1, tuple(negative): -1}


def conic_monomial(left, right, variables):
    answer = {(0,) * variables: 1}
    for triple in (left, right):
        for i, j in combinations(triple, 2):
            answer = multiply_polynomials(answer, bracket(i, j, variables))
    return answer


def check_symbolic_relations():
    r1_monomials = [
        ((0, 1, 2), (0, 3, 4)),
        ((0, 1, 3), (0, 2, 4)),
        ((0, 1, 4), (0, 2, 3)),
    ]
    r1_polynomials = [conic_monomial(*pair, 5) for pair in r1_monomials]
    relation = {}
    for coefficient, polynomial in zip((1, -1, 1), r1_polynomials):
        relation = add_polynomials(relation, polynomial, coefficient)
    assert relation == {}

    r0_monomials = []
    for left in combinations(range(6), 3):
        right = tuple(i for i in range(6) if i not in left)
        if left < right:
            r0_monomials.append((left, right))
    r0_polynomials = [conic_monomial(*pair, 6) for pair in r0_monomials]
    relations = [
        (-1, 1, -1, 1, 0, 0, 0, 0, 0, 0),
        (1, 0, 0, 0, 1, -1, 1, 0, 0, 0),
        (1, -1, 1, 0, 1, -1, 0, 1, 0, 0),
        (1, 0, 1, 0, 0, -1, 0, 0, 1, 0),
        (-1, 1, 0, 0, -1, 0, 0, 0, 0, 1),
    ]
    for coefficients in relations:
        relation = {}
        for coefficient, polynomial in zip(coefficients, r0_polynomials):
            relation = add_polynomials(relation, polynomial, coefficient)
        assert relation == {}
    return r1_monomials, r0_monomials


def check_boundary_witnesses(r1_monomials, r0_monomials):
    standard = tuple(range(1, 9))
    r1_rows = [
        primitive_row(evaluations(r1_monomials, {0, 1, 2, 5}, standard)),
        primitive_row(evaluations(r1_monomials, {0, 1, 3, 5}, standard)),
    ]
    assert r1_rows == [(0, 1, 1), (1, 0, -1)]
    assert determinant([[row[j] for j in (0, 1)] for row in r1_rows]) == -1

    specifications = [
        ({0, 1, 2, 3}, standard),
        ({0, 1, 2, 4}, standard),
        ({0, 1, 2, 5}, standard),
        ({0, 1, 3, 4}, standard),
        ({0, 1, 2, 3}, (5, 2, 3, 4, 5, 6, 7, 8)),
    ]
    r0_rows = [
        primitive_row(evaluations(r0_monomials, cut, parameters))
        for cut, parameters in specifications
    ]
    expected = [
        (0, 0, 1, 1, 0, 4, 4, 3, 3, 0),
        (0, 1, 0, -1, 3, 0, -3, -2, 0, 2),
        (0, -3, -3, 0, -8, -8, 0, 0, -5, -5),
        (1, 0, 0, 1, -9, -8, 0, 0, -9, -8),
        (0, 0, -3, -3, 0, -4, -4, -1, -1, 0),
    ]
    assert r0_rows == expected
    columns = (0, 1, 2, 4, 5)
    assert determinant([[row[j] for j in columns] for row in r0_rows]) == 8


def check_weight_enumeration():
    groups = {}
    for position, left in enumerate(TRIPLES):
        for right in TRIPLES[position:]:
            weight = tuple(int(i in left) + int(i in right) for i in range(8))
            groups.setdefault(weight, []).append((left, right))
    histogram = {}
    for monomials in groups.values():
        histogram[len(monomials)] = histogram.get(len(monomials), 0) + 1
    assert len(groups) == 784
    assert histogram == {1: 476, 3: 280, 10: 28}


def main():
    r1_monomials, r0_monomials = check_symbolic_relations()
    check_boundary_witnesses(r1_monomials, r0_monomials)
    check_weight_enumeration()
    print("PASS: symbolic conic identities and exact balanced-boundary ranks agree.")


if __name__ == "__main__":
    main()
