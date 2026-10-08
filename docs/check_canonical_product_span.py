"""Exact span test for canonical weighted source products.

The example deliberately has a common factor 3 in all source determinants.
It gives a nonsaturated product span inside the full contact-jet lattice.
"""

from itertools import combinations, combinations_with_replacement
from math import gcd
from fractions import Fraction

from check_higher_degree_gram_jet_barrier import (
    add_congruence, evaluate, transform_matrix,
)
from check_oriented_contact_lattice_rigidity import norm, mul, conj


def polynomial_product(rows, multiplicities):
    coeffs = [1]
    for (r, s), multiplicity in zip(rows, multiplicities):
        linear = [-s, r]  # det((r,s),(X,Y))
        for _ in range(multiplicity):
            out = [0] * (len(coeffs) + 1)
            for i, a in enumerate(coeffs):
                out[i] += a * linear[0]
                out[i + 1] += a * linear[1]
            coeffs = out
    return coeffs


def determinant(matrix):
    matrix = [row[:] for row in matrix]
    n = len(matrix)
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if matrix[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
            sign = -sign
        value = matrix[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = matrix[i][j] * value - matrix[i][k] * matrix[k][j]
                assert numerator % previous == 0
                matrix[i][j] = numerator // previous
        previous = value
    return sign * matrix[-1][-1]


def gcd_of_maximal_minors(columns):
    size = len(columns[0])
    value = 0
    for indices in combinations(range(len(columns)), size):
        matrix = [[columns[j][i] for j in indices] for i in range(size)]
        value = gcd(value, abs(determinant(matrix)))
    return value


def solve_columns(basis, vector):
    """Exact solve basis * x = vector for the small square matrices here."""
    n = len(basis)
    matrix = [[Fraction(basis[i][j]) for j in range(n)] + [Fraction(vector[i])]
              for i in range(n)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if matrix[i][col])
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        divisor = matrix[col][col]
        matrix[col] = [x / divisor for x in matrix[col]]
        for i in range(n):
            if i == col:
                continue
            factor = matrix[i][col]
            matrix[i] = [x - factor * y for x, y in zip(matrix[i], matrix[col])]
    return [matrix[i][-1] for i in range(n)]


def gamma_basis(degree, jet, contacts):
    basis = [[int(i == j) for j in range(degree + 1)]
             for i in range(degree + 1)]
    for modulus, lam in contacts:
        matrix = transform_matrix(degree, lam)
        for k in range(jet):
            add_congruence(basis, matrix[k], modulus ** (jet - k))
    return [[basis[i][j] for j in range(degree + 1)]
            for i in range(degree + 1)]


def gaussian_gcd(values):
    def rounded_quotient(a, b):
        z = mul(a, conj(b))
        denominator = norm(b)
        return ((2 * z[0] + denominator) // (2 * denominator),
                (2 * z[1] + denominator) // (2 * denominator))

    result = (0, 0)
    for value in values:
        a, b = value, result
        while b != (0, 0):
            q = rounded_quotient(a, b)
            qb = mul(q, b)
            a, b = b, (a[0] - qb[0], a[1] - qb[1])
        result = a
    return result


def check_example():
    # U is the identity frame, so Q_U=X^2+Y^2.  Each contact generator is
    # also its source row; its odd split norm is pairwise coprime.
    rows = [(1, 2), (4, 5), (7, 8)]
    contacts = [(5, (1, 2)), (41, (4, 5)), (113, (7, 8))]
    assert all(gcd(*row) == 1 for row in rows)
    assert [r * s2 - s * r2 for (r, s), (r2, s2)
            in combinations(rows, 2)] == [-3, -6, -3]
    D = 5 * 41 * 113

    for jet in (1, 2):
        for degree in (1, 2, 3):
            products, values = [], []
            for cuts in combinations_with_replacement(range(degree + 1), 2):
                multiplicities = (cuts[0], cuts[1] - cuts[0], degree - cuts[1])
                multiplier = 1
                for index, (modulus, _) in enumerate(contacts):
                    multiplier *= modulus ** max(jet - multiplicities[index], 0)
                polynomial = polynomial_product(rows, multiplicities)
                polynomial = [multiplier * x for x in polynomial]
                polynomial += [0] * (degree + 1 - len(polynomial))
                products.append(polynomial)
                values.append(evaluate(polynomial, (1, 0), (0, 1)))

            full_basis = gamma_basis(degree, jet, contacts)
            coordinates = []
            for polynomial in products:
                solution = solve_columns(full_basis, polynomial)
                assert all(x.denominator == 1 for x in solution)
                coordinates.append([int(x) for x in solution])

            full_index = D ** (jet * (jet + 1) // 2)
            product_index = gcd_of_maximal_minors(products)
            relative_index = gcd_of_maximal_minors(coordinates)
            assert product_index == full_index * relative_index
            assert relative_index == 3 ** (degree * (degree + 1) // 2)
            evaluation_norm = norm(gaussian_gcd(values))
            expected_evaluation_norm = D if jet == 1 else D ** min(2, degree + 2)
            # For (h,d)=(2,1), the canonical image has one extra D factor;
            # for d>=2 it reaches the expected image ideal (Delta^2).
            if jet == 2 and degree == 1:
                expected_evaluation_norm = D ** 3
            assert evaluation_norm == expected_evaluation_norm
            print(f"jet={jet}, degree={degree}: [Gamma:Lambda]={relative_index}, eval-norm={evaluation_norm}")


if __name__ == "__main__":
    check_example()
    print("Canonical products are nonsaturated; their Gaussian image can be nonsaturated as well.")
