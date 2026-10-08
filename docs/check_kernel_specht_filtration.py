"""Finite exact checks for the smaller-cut Specht filtration.

Checks Gram-minor null-collapse degrees for small highest weights,
Specht hook dimensions, termination, and the successive-minimum index
arithmetic. The general filtration uses polarization and representation
theory, not these finite examples.
"""

from itertools import permutations
from fractions import Fraction
from math import comb, factorial


ZERO = (0, 0, 0)


def add(*polynomials):
    out = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            out[exponent] = out.get(exponent, 0) + coefficient
    return {exponent: value for exponent, value in out.items() if value}


def scale(polynomial, scalar):
    return {exponent: scalar * value for exponent, value in polynomial.items()
            if scalar * value}


def multiply(left, right):
    out = {}
    for alpha, a in left.items():
        for beta, b in right.items():
            gamma = tuple(alpha[i] + beta[i] for i in range(3))
            out[gamma] = out.get(gamma, 0) + a * b
    return {exponent: value for exponent, value in out.items() if value}


def power(polynomial, exponent):
    result = {ZERO: 1}
    for _ in range(exponent):
        result = multiply(result, polynomial)
    return result


def degree(polynomial):
    return max(sum(exponent) for exponent in polynomial)


def top(polynomial):
    d = degree(polynomial)
    return {exponent: coefficient for exponent, coefficient in polynomial.items()
            if sum(exponent) == d}


def dot(x, y):
    return x[0] * y[2] + x[2] * y[0] - 2 * x[1] * y[1]


def determinant(matrix):
    size = len(matrix)
    out = {}
    for ordering in permutations(range(size)):
        inversions = sum(ordering[i] > ordering[j]
                         for i in range(size) for j in range(i + 1, size))
        term = {ZERO: (-1) ** inversions}
        for i, j in enumerate(ordering):
            term = multiply(term, matrix[i][j])
        out = add(out, term)
    return out


def scalar_determinant(columns):
    a, b, c = columns
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - b[0] * (a[1] * c[2] - a[2] * c[1])
            + c[0] * (a[1] * b[2] - a[2] * b[1]))


def gram_minors():
    u = (1, 0, 0)
    assert dot(u, u) == 0
    w = ((0, 0, 1), (0, 1, 1), (1, 1, 1))
    p = tuple(dot(u, vector) for vector in w)
    variables = [{tuple(int(i == j) for i in range(3)): 1}
                 for j in range(3)]
    gram = [[add({ZERO: dot(w[i], w[j])},
                 scale(variables[i], p[j]), scale(variables[j], p[i]))
             for j in range(3)] for i in range(3)]
    minors = [determinant([row[:size] for row in gram[:size]])
              for size in (1, 2, 3)]

    expected1 = scale(variables[0], 2 * p[0])
    expected2 = scale(power(add(scale(variables[0], p[1]),
                                scale(variables[1], -p[0])), 2), -1)
    cofactors = []
    for i in range(3):
        columns = list(w)
        columns[i] = u
        cofactors.append(scalar_determinant(columns))
    leading_det = add(*(scale(variables[i], cofactors[i]) for i in range(3)))
    bilinear_matrix = ((dot((1, 0, 0), (1, 0, 0)),
                        dot((1, 0, 0), (0, 1, 0)),
                        dot((1, 0, 0), (0, 0, 1))),
                       (dot((0, 1, 0), (1, 0, 0)),
                        dot((0, 1, 0), (0, 1, 0)),
                        dot((0, 1, 0), (0, 0, 1))),
                       (dot((0, 0, 1), (1, 0, 0)),
                        dot((0, 0, 1), (0, 1, 0)),
                        dot((0, 0, 1), (0, 0, 1))))
    det_form = scalar_determinant(bilinear_matrix)
    expected3 = scale(power(leading_det, 2), det_form)
    assert [degree(minor) for minor in minors] == [1, 2, 2]
    assert [top(minor) for minor in minors] == [expected1, expected2, expected3]
    return minors


def triples(q):
    for a in range(q, -1, -1):
        for b in range(a, -1, -1):
            c = q - a - b
            if 0 <= c <= b:
                yield a, b, c


def specht_dimension(a, b, c):
    m = 2 * (a + b + c)
    numerator = (factorial(m) * (2 * a - 2 * b + 1)
                 * (2 * a - 2 * c + 2) * (2 * b - 2 * c + 1))
    denominator = factorial(2 * a + 2) * factorial(2 * b + 1) * factorial(2 * c)
    assert numerator % denominator == 0
    return numerator // denominator


def invariant_dimension_by_weights(m):
    coefficients = [1]
    for _ in range(m):
        next_coefficients = [0] * (len(coefficients) + 2)
        for exponent, coefficient in enumerate(coefficients):
            for shift in range(3):
                next_coefficients[exponent + shift] += coefficient
        coefficients = next_coefficients
    return coefficients[m] - coefficients[m - 1]


def evaluation_exponent(m, k):
    q = m // 2
    base = m * 2 ** (m - 2) - q * comb(m, q)
    discount = sum(comb(m, s) * max(0, k - abs(s - q))
                   for s in range(m + 1))
    return base - discount


def profile(m):
    q = m // 2
    dimensions = {k: sum(specht_dimension(a, b, c)
                         for a, b, c in triples(q) if c >= k)
                  for k in range(q + 2)}
    assert dimensions[0] == invariant_dimension_by_weights(m)
    assert dimensions[q // 3 + 1] == 0
    assert dimensions[q // 3] > 0
    for k in range(q // 3 + 1):
        layer = sum(specht_dimension(a, b, c)
                    for a, b, c in triples(q) if c == k)
        assert dimensions[k] - dimensions[k + 1] == layer
        assert layer > 0
        if k >= 1 and dimensions[k] > 1:
            rank = dimensions[k] - 1
            next_rank = max(0, dimensions[k + 1] - 1)
            if dimensions[k + 1]:
                first_forced_index = next_rank + 1
                denominator = rank - first_forced_index + 1
                assert first_forced_index == dimensions[k + 1]
                assert denominator == layer
            else:
                assert rank == dimensions[k] - 1
    return dimensions


def main():
    minors = gram_minors()
    checked_weights = 0
    for q in range(1, 9):
        for a, b, c in triples(q):
            highest = multiply(multiply(power(minors[0], a - b),
                                        power(minors[1], b - c)),
                               power(minors[2], c))
            assert degree(highest) == a + b
            checked_weights += 1

    dimensions = {m: profile(m) for m in range(2, 22, 2)}
    assert dimensions[6][1] == 5
    assert dimensions[8][1] == 56
    assert dimensions[10][1] == 477 and dimensions[10][2] == 0
    assert dimensions[12][2] == 462
    known_exponents = {(6, 1): 16, (8, 1): 162, (10, 1): 1048,
                       (12, 1): 5820, (12, 2): 3312, (18, 3): 357528}
    for key, expected in known_exponents.items():
        assert evaluation_exponent(*key) == expected
    ratios = {}
    for m, ranks in dimensions.items():
        for k in range(1, m // 6 + 1):
            rank, next_rank = ranks[k], ranks[k + 1]
            if rank <= 1:
                continue
            denominator = rank - next_rank if next_rank else rank - 1
            assert denominator > 0
            ratios[m, k] = Fraction(evaluation_exponent(m, k), denominator)
    assert ratios[12, 1] == Fraction(5820, 3289)
    assert ratios[12, 2] == Fraction(3312, 461)
    print(f"PASS: Gram-minor top degrees for {checked_weights} highest weights;"
          " SO(3) Specht dimensions agree with spin-one weight counts.")
    print("      K1 dimensions m=6,8,10: 5,56,477; K2(10)=0, K2(12)=462.")
    print("      Layer sizes, termination, height coefficients, and minima indices checked through m=20.")


if __name__ == "__main__":
    main()
