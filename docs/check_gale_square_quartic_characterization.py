#!/usr/bin/env python3
"""Exact checks for Gale interpolation and the square-quartic obstruction."""

from fractions import Fraction
from functools import reduce
from itertools import combinations, combinations_with_replacement
from math import gcd, isqrt, lcm


ZERO = Fraction(0)
ONE = Fraction(1)
ELLIPTIC_POINT = (Fraction(3), Fraction(5))


def elliptic_add(left, right):
    """Add points on y^2=x^3-2; None is the point at infinity."""
    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right and y_left == -y_right:
        return None
    if left == right:
        if y_left == 0:
            return None
        slope = 3 * x_left * x_left / (2 * y_left)
    else:
        slope = (y_right - y_left) / (x_right - x_left)
    x_sum = slope * slope - x_left - x_right
    y_sum = slope * (x_left - x_sum) - y_left
    return x_sum, y_sum


def first_seven_multiples():
    points = []
    point = None
    for _ in range(7):
        point = elliptic_add(point, ELLIPTIC_POINT)
        assert point is not None
        x, y = point
        assert y * y == x * x * x - 2
        points.append(point)
    assert len({x for x, _ in points}) == 7
    return points


def polynomial_add(left, right):
    size = max(len(left), len(right))
    answer = [ZERO] * size
    for i in range(size):
        answer[i] = (left[i] if i < len(left) else ZERO) + (
            right[i] if i < len(right) else ZERO
        )
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def polynomial_scale(polynomial, scalar):
    return [scalar * coefficient for coefficient in polynomial]


def polynomial_multiply(left, right):
    answer = [ZERO] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            answer[i + j] += left_coefficient * right_coefficient
    return answer


def interpolate(nodes, values):
    answer = [ZERO]
    for i, node in enumerate(nodes):
        numerator = [ONE]
        denominator = ONE
        for j, other in enumerate(nodes):
            if i != j:
                numerator = polynomial_multiply(numerator, [-other, ONE])
                denominator *= node - other
        answer = polynomial_add(answer, polynomial_scale(numerator, values[i] / denominator))
    return answer


def primitive_pluckers(rows):
    rational = {}
    for i, j in combinations(range(7), 2):
        rational[i, j] = rows[i][0] * rows[j][1] - rows[i][1] * rows[j][0]
    common_denominator = 1
    for value in rational.values():
        common_denominator = lcm(common_denominator, value.denominator)
    integral = {pair: int(value * common_denominator) for pair, value in rational.items()}
    content = 0
    for value in integral.values():
        content = gcd(content, abs(value))
    primitive = {pair: value // content for pair, value in integral.items()}
    assert reduce(gcd, (abs(value) for value in primitive.values())) == 1
    return rational, primitive


def signed_coordinate(pluckers, i, j):
    return pluckers[i, j] if i < j else -pluckers[j, i]


def is_rational_square(value):
    if value < 0:
        return False
    return (
        isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


def check_gale_example():
    points = first_seven_multiples()
    nodes = [x for x, _ in points]
    ordinates = [y for _, y in points]
    derivatives = [
        reduce(
            lambda left, right: left * right,
            (nodes[i] - nodes[j] for j in range(7) if i != j),
            ONE,
        )
        for i in range(7)
    ]
    values = [y * y for y in ordinates]
    scalings = [values[i] / derivatives[i] for i in range(7)]
    rows = [(scalings[i], scalings[i] * nodes[i]) for i in range(7)]
    assert sum(row[0] for row in rows) == 0
    assert sum(row[1] for row in rows) == 0

    interpolant = interpolate(nodes, [scalings[i] * derivatives[i] for i in range(7)])
    assert interpolant == [Fraction(-2), ZERO, ZERO, ONE]
    assert (len(interpolant) - 1) % 2 == 1

    rational, primitive = primitive_pluckers(rows)
    for i, j, k, ell in combinations(range(7), 4):
        assert (
            primitive[i, j] * primitive[k, ell]
            - primitive[i, k] * primitive[j, ell]
            + primitive[i, ell] * primitive[j, k]
            == 0
        )

    product_scalings = reduce(lambda left, right: left * right, scalings, ONE)
    for i in range(7):
        rational_star = reduce(
            lambda left, j: left * signed_coordinate(rational, i, j),
            (j for j in range(7) if j != i),
            ONE,
        )
        expected = product_scalings * values[i] ** 5 / derivatives[i] ** 4
        assert rational_star == expected
        assert is_rational_square(-rational_star)

        primitive_star = reduce(
            lambda left, j: left * signed_coordinate(primitive, i, j),
            (j for j in range(7) if j != i),
            1,
        )
        assert primitive_star < 0
        assert isqrt(-primitive_star) ** 2 == -primitive_star

    return primitive


def monomial_multiply(left, right):
    answer = {}
    for left_powers, left_coefficient in left.items():
        for right_powers, right_coefficient in right.items():
            powers = tuple(a + b for a, b in zip(left_powers, right_powers))
            answer[powers] = answer.get(powers, 0) + left_coefficient * right_coefficient
    return answer


def determinant(matrix):
    matrix = [list(map(Fraction, row)) for row in matrix]
    answer = ONE
    for column in range(len(matrix)):
        pivot = next(
            row
            for row in range(column, len(matrix))
            if matrix[row][column]
        )
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            answer = -answer
        pivot_value = matrix[column][column]
        answer *= pivot_value
        for row in range(column + 1, len(matrix)):
            multiplier = matrix[row][column] / pivot_value
            for j in range(column, len(matrix)):
                matrix[row][j] -= multiplier * matrix[column][j]
    return answer


def check_square_quartic_quadratic_rank():
    # Variables are A,B,C.  These are the five descending coefficients.
    coefficients = [
        {(2, 0, 0): 1},
        {(1, 1, 0): 2},
        {(0, 2, 0): 1, (1, 0, 1): 2},
        {(0, 1, 1): 2},
        {(0, 0, 2): 1},
    ]
    products = [
        monomial_multiply(coefficients[i], coefficients[j])
        for i, j in combinations_with_replacement(range(5), 2)
    ]
    degree_four_monomials = [
        (4 - exponent_b - exponent_c, exponent_b, exponent_c)
        for exponent_b in range(5)
        for exponent_c in range(5 - exponent_b)
    ]
    coefficient_matrix = [
        [product.get(monomial, 0) for product in products]
        for monomial in degree_four_monomials
    ]
    assert determinant(coefficient_matrix) == -16384


def main():
    primitive = check_gale_example()
    check_square_quartic_quadratic_rank()
    largest = max(len(str(abs(value))) for value in primitive.values())
    print(
        "PASS: seven exact rows have primitive negative-square stars, "
        f"a nonsquare cubic interpolant, and Pluckers of at most {largest} digits."
    )
    print("PASS: the square-quartic coefficient map has quadratic rank 15.")


if __name__ == "__main__":
    main()
