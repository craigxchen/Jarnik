"""Exact Gram compression, weighted contents, and normalized constants."""

from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd, lcm


def prod(values):
    return reduce(lambda a, b: a * b, values, F(1))


def primitive(values):
    denominator = lcm(*(x.denominator for x in values))
    integer = [int(x * denominator) for x in values]
    common = reduce(gcd, map(abs, integer))
    return [x // common for x in integer]


def inverse_and_determinant(matrix):
    size = len(matrix)
    work = [[F(x) for x in row] + [F(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    determinant = F(1)
    for j in range(size):
        pivot = next(i for i in range(j, size) if work[i][j])
        if pivot != j:
            work[j], work[pivot] = work[pivot], work[j]
            determinant = -determinant
        factor = work[j][j]
        determinant *= factor
        work[j] = [x / factor for x in work[j]]
        for i in range(size):
            if i != j:
                factor = work[i][j]
                work[i] = [x - factor * y for x, y in zip(work[i], work[j])]
    return [row[size:] for row in work], determinant


def weighted_matrix(nodes):
    differences = [prod(t - r for j, r in enumerate(nodes) if j != i)
                   for i, t in enumerate(nodes)]
    return [[t ** k / difference for k in range(5)]
            for t, difference in zip(nodes, differences)], differences


def clear_matrix(matrix):
    scale = lcm(*(x.denominator for row in matrix for x in row))
    return scale, [[int(x * scale) for x in row] for row in matrix]


def check_base(nodes):
    nodes = list(map(F, nodes))
    evaluation, differences = weighted_matrix(nodes)
    scale, integer_matrix = clear_matrix(evaluation)
    inverse, determinant = inverse_and_determinant(integer_matrix[:5])
    determinant = abs(int(determinant))
    inverse_norm = max(sum(map(abs, row)) for row in inverse)
    constant = 1 / (4 * determinant * inverse_norm)
    count = positive_count = 0
    for a, b, c, d in product(range(-2, 3), repeat=4):
        if a * d == b * c or any(c * t + d == 0 for t in nodes):
            continue
        A, B, C = a * a + c * c, a * b + c * d, b * b + d * d
        common = gcd(gcd(A, B), C)
        A, B, C = A // common, B // common, C // common
        H = max(abs(A), abs(B), abs(C))
        J = A * C - B * B
        assert J > 0 and (a * d - b * c) ** 2 == J * common ** 2
        coefficients = [C * C, 4 * B * C, 2 * A * C + 4 * B * B,
                        4 * A * B, A * A]
        assert reduce(gcd, map(abs, coefficients)) <= 4
        assert max(map(abs, coefficients)) >= H * H
        raw_integer = [sum(x * y for x, y in zip(row, coefficients))
                       for row in integer_matrix]
        content = reduce(gcd, map(abs, raw_integer))
        assert (4 * determinant) % content == 0
        actual = [(a * t + b) / (c * t + d) for t in nodes]
        _, actual_differences = weighted_matrix(actual)
        actual_weights = primitive([(t * t + 1) ** 2 / diff
                                    for t, diff in zip(actual, actual_differences)])
        gram_weights = primitive([F(x) for x in raw_integer])
        assert actual_weights == gram_weights or actual_weights == [-x for x in gram_weights]
        U = max(map(abs, gram_weights))
        assert U >= constant * H * H
        if min(actual) > 0:
            X = min(actual)
            i, j = 0, 1
            base_chord = 4 * (nodes[i] - nodes[j]) ** 2 / (
                (nodes[i] ** 2 + 1) * (nodes[j] ** 2 + 1))
            image_chord = 4 * (actual[i] - actual[j]) ** 2 / (
                (actual[i] ** 2 + 1) * (actual[j] ** 2 + 1))
            assert image_chord >= J * base_chord / (A + C) ** 2
            assert image_chord <= 4 / (X * X)
            assert U >= constant * base_chord * J * X * X / 16
            positive_count += 1
        count += 1
    return count, positive_count


def normalized_constant(nodes, expected_scale, expected_constant):
    p, q, r = map(F, nodes[:3])
    finite = [F(0), F(1)] + [(F(t) - p) * (q - r) /
                              ((F(t) - r) * (q - p)) for t in nodes[3:]]
    evaluation, _ = weighted_matrix(finite)
    scale, matrix = clear_matrix(evaluation)
    inverse, determinant = inverse_and_determinant(matrix)
    inverse_norm = max(sum(map(abs, row)) for row in inverse)
    constant = 1 / (4 * abs(determinant) * inverse_norm)
    assert scale == expected_scale
    assert constant == expected_constant
    # The infinity weight is minus the leading coefficient of the quartic.
    for coefficients in ([1, 2, 3, 4, 5], [4, 12, 13, 6, 1]):
        finite_values = [sum(x * y for x, y in zip(row, coefficients))
                         for row in evaluation]
        assert sum(finite_values) == coefficients[-1]


def main():
    counts = [check_base(nodes) for nodes in
              ([2, 3, 5, 12, 17, 23], [1, 2, 3, 4, 7, 11],
               [F(1, 2), F(2, 3), F(4, 5), F(6, 7), F(8, 9), F(10, 11)])]
    normalized_constant([1, 2, 3, 4, 7, 11], 315, F(1, 9601804800))
    normalized_constant([2, 3, 5, 12, 17, 23], 4950, F(1, 4034503242000000))
    print("PASS:", sum(n for n, _ in counts), "rational Mobius images;",
          sum(n for _, n in counts), "positive-image chord/height inequalities")
    print("PASS: two exact normalized evaluation constants and infinity rows")


if __name__ == "__main__":
    main()
