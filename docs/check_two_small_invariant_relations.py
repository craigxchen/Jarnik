"""Exact certificate for the two-small-relations theorem on six rows.

The calculation concerns polynomial restrictions, not numerical sample
points.  The arithmetic application is proved in two_small_invariant_relations.md.
Only the Python standard library is used.
"""

from fractions import Fraction
from itertools import combinations

from check_near_balanced_invariant_congruences import (
    ZERO, add, bracket, mul, quadratic_coefficient, scale, variable,
)


TRIPLES = ((0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 2, 3), (0, 2, 4))
PAIRS = tuple(combinations(range(6), 2))
WEDGE_COLUMNS = tuple(combinations(range(5), 2))
MINOR_LABELS = ("12", "13", "14", "15", "23", "24", "25", "34", "35", "45")


def triangle_restriction(triple, inside):
    result = {ZERO: 1}
    complement = tuple(j for j in range(6) if j not in triple)
    for a, b, c in (triple, complement):
        for i, j in ((a, b), (b, c), (c, a)):
            result = mul(result, bracket(i, j, inside))
    return result


def determinant(matrix):
    rows = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next((j for j in range(column, len(rows)) if rows[j][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[pivot], rows[column] = rows[column], rows[pivot]
            result = -result
        value = rows[column][column]
        result *= value
        for j in range(column + 1, len(rows)):
            factor = rows[j][column] / value
            rows[j] = [a - factor * b for a, b in zip(rows[j], rows[column])]
    return result


def main():
    for inside in combinations(range(6), 3):
        for triple in TRIPLES:
            assert not triangle_restriction(triple, set(inside))

    restriction_rows = {}
    wedge_rows = {}
    for pair in PAIRS:
        inside = set(pair)
        a, b, c, d = (j for j in range(6) if j not in inside)
        first = mul(add(variable(a), scale(variable(b), -1)),
                    add(variable(c), scale(variable(d), -1)))
        second = mul(add(variable(a), scale(variable(c), -1)),
                     add(variable(b), scale(variable(d), -1)))
        coefficients = []
        for triple in TRIPLES:
            polynomial = triangle_restriction(triple, inside)
            alpha = quadratic_coefficient(polynomial, a, c)
            beta = quadratic_coefficient(polynomial, a, b)
            assert polynomial == add(scale(first, alpha), scale(second, beta))
            coefficients.append((alpha, beta))
        label = "".join(str(i + 1) for i in pair)
        restriction_rows[label] = coefficients
        wedge_rows[label] = [
            coefficients[i][0] * coefficients[j][1]
            - coefficients[i][1] * coefficients[j][0]
            for i, j in WEDGE_COLUMNS
        ]

    minor = [wedge_rows[label] for label in MINOR_LABELS]
    assert determinant(minor) == -2
    print("All 5 triangle products vanish at every balanced restriction.")
    print("All 15 first-smaller restrictions match their two matching quadratics exactly.")
    print("Selected rows:", ", ".join(MINOR_LABELS))
    print("Their 10 by 10 wedge minor has determinant -2; wedge map rank = 10.")
    print("Restriction coefficients (five pairs in the documented basis):")
    for label, coefficients in restriction_rows.items():
        print(label, coefficients)
    print("No numerical zero relation or lattice construction is asserted by this checker.")


if __name__ == "__main__":
    main()
