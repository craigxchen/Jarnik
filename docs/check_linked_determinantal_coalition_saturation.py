"""Small exact audits of the linked-ideal and invariant-lattice arguments.

The all-residue dimension and flatness proof is in the companion note.
This checker adds no new large source or twelve-row rank fixture.
"""

from fractions import Fraction
from itertools import combinations
from math import lcm
from random import Random

from check_coherent_kernel_eagon_northcott_resolution import det_mod
from check_distinct_configuration_determinantal_kernel_generation import source_basis


def inverse(matrix):
    size = len(matrix)
    work = [[Fraction(value) for value in row]
            + [Fraction(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(i for i in range(column, size) if work[i][column])
        work[pivot], work[column] = work[column], work[pivot]
        divisor = work[column][column]
        work[column] = [value / divisor for value in work[column]]
        for i in range(size):
            if i != column:
                multiplier = work[i][column]
                work[i] = [left - multiplier * right
                           for left, right in zip(work[i], work[column])]
    return [row[size:] for row in work]


def determinant_integer(matrix):
    matrix = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(matrix)):
        pivot = next(i for i in range(column, len(matrix)) if matrix[i][column])
        if pivot != column:
            matrix[pivot], matrix[column] = matrix[column], matrix[pivot]
            result *= -1
        value = matrix[column][column]
        result *= value
        for i in range(column + 1, len(matrix)):
            multiplier = matrix[i][column] / value
            for j in range(column + 1, len(matrix)):
                matrix[i][j] -= multiplier * matrix[column][j]
    return result


def check_invariant_lattice():
    basis = source_basis(3, 2, tuple(range(6)))
    assert len(basis) == 5
    leading = [min(polynomial) for polynomial in basis]
    assert len(set(leading)) == 5
    minor = [[basis[j].get(leading[i], 0) for j in range(5)] for i in range(5)]
    assert abs(determinant_integer(minor)) == 1

    # The ordinary tensor inner product is SU(3)-invariant. Thus this Gram
    # formula gives the same invariant projection as the Reynolds operator.
    gram = [[sum(value * basis[j].get(key, 0)
                 for key, value in basis[i].items())
             for j in range(5)] for i in range(5)]
    gram_inverse = inverse(gram)
    keys = sorted(set().union(*(polynomial.keys() for polynomial in basis)))
    denominator = 1
    for left in keys:
        for right in keys:
            coefficient = sum(basis[i].get(left, 0) * gram_inverse[i][j]
                              * basis[j].get(right, 0)
                              for i in range(5) for j in range(5))
            denominator = lcm(denominator, coefficient.denominator)
    assert denominator == 72
    assert 144 % denominator == 0

    # Existing G_A,...,G_E relative to T_123,...,T_135.
    change = ((2, -1, 1, 1, -1), (1, -1, 0, 1, 0),
              (1, -1, 1, 0, 0), (1, -1, 0, 0, 0), (1, 0, 0, 0, 0))
    assert abs(determinant_integer(change)) == 1
    return denominator


def interpolate_mod(values, prime):
    """Values at 0,...,n determine the coefficients of degree at most n."""
    degree = len(values) - 1
    result = [0] * (degree + 1)
    for node, value in enumerate(values):
        polynomial = [1]
        denominator = 1
        for other in range(degree + 1):
            if other == node:
                continue
            expanded = [0] * (len(polynomial) + 1)
            for power, coefficient in enumerate(polynomial):
                expanded[power] -= other * coefficient
                expanded[power + 1] += coefficient
            polynomial = expanded
            denominator *= node - other
        scale = value * pow(denominator % prime, -1, prime)
        for power, coefficient in enumerate(polynomial):
            result[power] = (result[power] + scale * coefficient) % prime
    return result


def check_linked_minors():
    prime = 101
    random = Random(20261003)
    tested = 0
    for n, m in ((2, 6), (3, 9)):
        for coalition_size in range(m + 1):
            vectors = [[random.randrange(prime) for _ in range(n)]
                       for _ in range(m)]
            slopes = [i + 1 for i in range(m)]

            def a_matrix(t, selected):
                return [vectors[i] + [entry * slopes[i]
                        * (t if i < coalition_size else 1) % prime
                        for entry in vectors[i]] for i in selected]

            def b_matrix(t, selected):
                return [[entry * (1 if i < coalition_size else t) % prime
                         for entry in vectors[i]]
                        + [entry * slopes[i] % prime for entry in vectors[i]]
                        for i in selected]

            nonzero_leading = 0
            for selected in combinations(range(m), 2 * n):
                hits = sum(i < coalition_size for i in selected)
                forced = max(0, hits - n)
                values = [det_mod(a_matrix(t, selected), prime)
                          for t in range(n + 1)]
                coefficients = interpolate_mod(values, prime)
                assert all(value == 0 for value in coefficients[:forced])
                specialized = det_mod(a_matrix(0, selected) if hits <= n
                                      else b_matrix(0, selected), prime)
                assert coefficients[forced] == specialized
                nonzero_leading += bool(specialized)
                for t in (2, 7):
                    assert det_mod(b_matrix(t, selected), prime) == (
                        pow(t, n - hits, prime) * det_mod(a_matrix(t, selected), prime)
                    ) % prime
                tested += 1
            assert nonzero_leading
    return tested


def main():
    denominator = check_invariant_lattice()
    print('PASS: determinant and gradient invariant lattices have unit minors.')
    print(f'PASS: six-row Reynolds denominator is {denominator}, dividing 144.')
    print(f'PASS: {check_linked_minors()} linked-minor and normalized-limit checks.')
    print('The uniform saturation theorem uses the linked dimension/flatness proof.')


if __name__ == '__main__':
    main()
