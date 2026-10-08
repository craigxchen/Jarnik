"""Finite exact audits for the all-distinct determinantal generation theorem.

The proof of primeness and all-configuration generation is in the companion
note. Finite modular ranks here are checks, not a substitute for that proof.
Only the Python standard library is required.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial

from check_local_tensor_coefficient_floor_sharpness import determinant_kernel, sign


PRIME = 101


def sparse_rank(vectors):
    pivots = {}
    for vector in vectors:
        vector = {key: int(value % PRIME) for key, value in vector.items()
                  if value % PRIME}
        while vector:
            lead = min(vector)
            if lead not in pivots:
                inverse = pow(vector[lead], -1, PRIME)
                pivots[lead] = {key: value * inverse % PRIME
                                for key, value in vector.items()}
                break
            factor = vector[lead]
            for key, value in pivots[lead].items():
                new = (vector.get(key, 0) - factor * value) % PRIME
                if new:
                    vector[key] = new
                else:
                    vector.pop(key, None)
    return len(pivots)


def matrix_rank(matrix):
    return sparse_rank(dict(enumerate(row)) for row in matrix)


def tableaux(n, d, rows=None):
    if rows is None:
        rows = ((),) * n
    label = sum(map(len, rows))
    if label == n * d:
        yield rows
        return
    for row in range(n):
        if len(rows[row]) == d or (row and len(rows[row - 1]) <= len(rows[row])):
            continue
        next_rows = list(rows)
        next_rows[row] += (label,)
        yield from tableaux(n, d, tuple(next_rows))


def determinant_products(n, labels, blocks):
    """Coefficients in row-color monomials, with rows ordered by labels."""
    positions = {label: i for i, label in enumerate(labels)}
    terms = {(-1,) * len(labels): 1}
    for block in blocks:
        next_terms = {}
        for colors, value in terms.items():
            for permutation in permutations(range(n)):
                target = list(colors)
                for label, color in zip(block, permutation):
                    target[positions[label]] = color
                next_terms[tuple(target)] = value * sign(permutation)
        terms = next_terms
    return terms


def source_basis(n, d, labels):
    if not d:
        return [{(): 1}]
    output = []
    for tableau in tableaux(n, d):
        blocks = [tuple(labels[i] for i in column) for column in zip(*tableau)]
        output.append(determinant_products(n, labels, blocks))
    return output


def coherent(polynomial, z):
    output = defaultdict(int)
    for colors, value in polynomial.items():
        key = []
        for i, color in enumerate(colors):
            if color == 1:
                value *= z[i]
            elif color >= 2:
                key.append((i, color))
        output[tuple(key)] += value
    return {key: value for key, value in output.items() if value}


def combine_disjoint(local, outside, selected, complement, m):
    answer = {}
    for local_colors, first in local.items():
        for outside_colors, second in outside.items():
            colors = [-1] * m
            for label, color in zip(selected, local_colors):
                colors[label] = color
            for label, color in zip(complement, outside_colors):
                colors[label] = color
            answer[tuple(colors)] = first * second
    return answer


def check_small_generation(n, d, z):
    m = n * d
    labels = tuple(range(m))
    basis = source_basis(n, d, labels)
    source_dimension = factorial(m)
    hook_product = 1
    for i in range(n):
        for j in range(d):
            hook_product *= n + d - i - j - 1
    source_dimension //= hook_product
    assert len(basis) == source_dimension
    assert sparse_rank(basis) == source_dimension
    evaluation_rank = sparse_rank(coherent(poly, z) for poly in basis)

    def generators():
        for selected in combinations(labels, 2 * n):
            complement = tuple(i for i in labels if i not in selected)
            local = determinant_kernel(n, [z[i] for i in selected])
            for outside in source_basis(n, d - 2, complement):
                poly = combine_disjoint(local, outside, selected, complement, m)
                assert not coherent(poly, z)
                yield poly

    generator_rank = sparse_rank(generators())
    assert generator_rank + evaluation_rank == source_dimension
    return source_dimension, evaluation_rank, generator_rank


def linear_combination(coefficients, polynomials):
    answer = defaultdict(int)
    for scalar, polynomial in zip(coefficients, polynomials):
        for monomial, value in polynomial.items():
            answer[monomial] += scalar * value
    return {key: value for key, value in answer.items() if value}


def check_six_row_polar(z):
    labels = tuple(range(6))
    triples = ((0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 2, 3), (0, 2, 4))
    pairs = [determinant_products(3, labels, (triple, tuple(i for i in labels
             if i not in triple))) for triple in triples]
    gradient_coefficients = ((2, -1, 1, 1, -1), (1, -1, 0, 1, 0),
                             (1, -1, 1, 0, 0), (1, -1, 0, 0, 0),
                             (1, 0, 0, 0, 0))
    gradients = [linear_combination(coefficients, pairs)
                 for coefficients in gradient_coefficients]
    matchings = (((0, 1), (2, 3), (4, 5)), ((0, 1), (2, 5), (3, 4)),
                 ((0, 3), (1, 2), (4, 5)), ((0, 5), (1, 2), (3, 4)),
                 ((0, 5), (1, 4), (2, 3)))
    matching_values = []
    for matching in matchings:
        value = 1
        for i, j in matching:
            value *= z[j] - z[i]
        matching_values.append(value)
    polar = linear_combination(matching_values, gradients)
    local = determinant_kernel(3, z)
    first = min(local)
    scalar = Fraction(polar[first], local[first])
    assert scalar
    assert polar == {key: scalar * value for key, value in local.items()}
    assert not coherent(polar, z)
    return scalar


def check_vandermonde_witness(n, m, z):
    exponents = (0, 1) + tuple(range(3, 2 * n - 2, 2))
    assert len(exponents) == n
    v = [[pow(t, exponent, PRIME) for exponent in exponents] for t in z]
    a = [row + [t * entry % PRIME for entry in row] for t, row in zip(z, v)]
    rank = matrix_rank(a)
    assert rank == 2 * n - 1
    # A*(e_2,-e_1)=0; each row normal is (-z_i,1,0,...,0).
    assert all(row[1] == row[n] for row in a)
    derivative = []
    for i, t in enumerate(z):
        row = [0] * (m * n)
        row[i * n], row[i * n + 1] = -t, 1
        derivative.append(row)
    assert matrix_rank(derivative) == m
    assert matrix_rank([x + y for x, y in zip(a, derivative)]) - rank == m - 2 * n + 1


def check_repeated_direction_scope():
    for n in range(2, 5):
        # Two groups of n directions: expected height, reducible determinant.
        z = (0,) * n + (1,) * n
        local = determinant_kernel(n, z)
        ordinary_product = determinant_products(
            n, tuple(range(2 * n)), (tuple(range(n)), tuple(range(n, 2 * n))))
        assert local == ordinary_product
        # One group of n+1: fewer than n outside rows, so the ideal is zero.
        assert not determinant_kernel(n, (0,) * (n + 1) + (1,) * (n - 1))


def main():
    witness_count = 0
    for n in range(2, 7):
        for d in range(2, 5):
            m = n * d
            check_vandermonde_witness(n, m, tuple(range(m)))
            check_vandermonde_witness(n, m, tuple(pow(i + 1, -1, PRIME)
                                                for i in range(m)))
            witness_count += 2
    print(f'{witness_count} Vandermonde rank and normal-map witnesses pass')

    fixtures = (tuple(range(6)), (-7, -2, 0, 3, 11, 29),
                tuple(Fraction(1, i + 1) for i in range(6)))
    scalars = [check_six_row_polar(z) for z in fixtures]
    print(f'six-row polar / local-determinant scalars: {scalars}')
    check_repeated_direction_scope()
    print('repeated-direction factorization and zero-ideal witnesses pass')

    for n, d in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2)):
        m = n * d
        for z in (tuple(range(m)), tuple(i * i + 3 for i in range(m))):
            ranks = check_small_generation(n, d, z)
            print(f'(n,d)=({n},{d}), z={z}: source/image/local ranks={ranks}')
    print('PASS: finite exact checks; the all-point theorem uses the scheme proof.')


if __name__ == '__main__':
    main()
