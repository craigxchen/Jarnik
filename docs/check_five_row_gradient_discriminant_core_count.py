"""Exact finite algebra for the five-row row-gradient core count.

The checker constructs literal bracket polynomials.  It verifies the
four-row graph expansion of each polarized row discriminant, then counts
the internal-edge minima at each Gaussian cut.  It makes no endpoint
height claim.
"""

from fractions import Fraction
from itertools import combinations_with_replacement
from math import lcm

from check_odd_support_relation_rigidity import graphs
from check_smaller_cut_relation_rank import add, derivative, mul, scale, variable


def bracket(n, i, j):
    xi, yi = variable(2 * n, 2 * i), variable(2 * n, 2 * i + 1)
    xj, yj = variable(2 * n, 2 * j), variable(2 * n, 2 * j + 1)
    return add(mul(xi, yj), scale(mul(yi, xj), -1))


def graph_poly(n, edges):
    poly = {(0,) * (2 * n): 1}
    for i, j in edges:
        poly = mul(poly, bracket(n, i, j))
    return poly


def reduce_with_pivots(poly, pivots):
    work = {m: Fraction(c) for m, c in poly.items()}
    for pivot, row in sorted(pivots.items()):
        factor = work.get(pivot, 0)
        if factor:
            work = add(work, scale(row, -factor))
    return work


def independent_graphs(n, degree):
    pivots = {}
    chosen = []
    for edges in graphs((degree,) * n):
        poly = graph_poly(n, edges)
        reduced = reduce_with_pivots(poly, pivots)
        if reduced:
            pivot = min(reduced)
            pivots[pivot] = scale(reduced, 1 / reduced[pivot])
            chosen.append((edges, poly))
    return chosen


def quadratic_coefficients(poly, i):
    parts = [{}, {}, {}]
    for monomial, coefficient in poly.items():
        xpower, ypower = monomial[2 * i : 2 * i + 2]
        assert xpower + ypower == 2
        reduced = monomial[:2 * i] + monomial[2 * i + 2 :]
        index = { (0, 2): 0, (1, 1): 1, (2, 0): 2 }[(xpower, ypower)]
        parts[index][reduced] = coefficient
    return parts  # c, b, a


def polarized_discriminant(first, second, same=False):
    c, b, a = first
    C, B, A = second
    if same:
        return add(mul(b, b), scale(mul(a, c), -4))
    return add(scale(mul(b, B), 2),
               scale(add(mul(a, C), mul(A, c)), -4))


def independent_monomials(basis):
    # Rows at independent monomials are found by taking successive ranks.
    monomials = sorted(set().union(*(poly.keys() for _, poly in basis)))
    chosen = []
    pivots = {}
    for monomial in monomials:
        vector = {j: poly.get(monomial, 0)
                  for j, (_, poly) in enumerate(basis) if poly.get(monomial, 0)}
        reduced = reduce_with_pivots(vector, pivots)
        if reduced:
            pivot = min(reduced)
            pivots[pivot] = scale(reduced, 1 / reduced[pivot])
            chosen.append(monomial)
        if len(chosen) == len(basis):
            break
    assert len(chosen) == len(basis)
    return chosen


def determinant(matrix):
    rows = [[Fraction(v) for v in row] for row in matrix]
    n = len(rows)
    det = Fraction(1)
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        if pivot != col:
            rows[col], rows[pivot] = rows[pivot], rows[col]
            det = -det
        factor = rows[col][col]
        det *= factor
        for i in range(col + 1, n):
            ratio = rows[i][col] / factor
            rows[i] = [x - ratio * y for x, y in zip(rows[i], rows[col])]
    return det


def inverse(matrix):
    rows = [[Fraction(v) for v in row] +
            [Fraction(int(i == j)) for j in range(len(matrix))]
            for i, row in enumerate(matrix)]
    n = len(rows)
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        factor = rows[col][col]
        rows[col] = [v / factor for v in rows[col]]
        for i in range(n):
            if i != col and rows[i][col]:
                factor = rows[i][col]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[col])]
    return [row[n:] for row in rows]


def solve_coefficients(basis, target):
    chosen = independent_monomials(basis)
    n = len(basis)
    rows = [[Fraction(poly.get(m, 0)) for _, poly in basis] +
            [Fraction(target.get(m, 0))] for m in chosen]
    for col in range(n):
        pivot = next(i for i in range(col, n) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        factor = rows[col][col]
        rows[col] = [v / factor for v in rows[col]]
        for i in range(n):
            if i != col and rows[i][col]:
                factor = rows[i][col]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[col])]
    coefficients = [row[-1] for row in rows]
    reconstructed = {}
    for coefficient, (_, poly) in zip(coefficients, basis):
        reconstructed = add(reconstructed, scale(poly, coefficient))
    assert reconstructed == target
    return coefficients


def internal_minimum(degree, a):
    inside = set(range(a))
    return min(sum(i in inside and j in inside for i, j in edges)
               for edges in graphs((degree,) * 4))


def incident_minima():
    output = {}
    for i_inside in (False, True):
        for a in range(5):
            inside = set(range(1, a + 1)) | ({0} if i_inside else set())
            values = []
            for edges in graphs((2,) * 5):
                internal = sum(i in inside and j in inside for i, j in edges)
                for i, j in edges:
                    if 0 not in (i, j):
                        continue
                    neighbor = j if i == 0 else i
                    removed = int(i_inside and neighbor in inside)
                    pi = internal - removed + int(neighbor in inside) - int(i_inside)
                    bar = internal - removed
                    values.append((pi, bar))
            output[(i_inside, a)] = (min(pi for pi, _ in values),
                                      min(bar for _, bar in values))
    return output


def evaluate(poly, values):
    result = 0
    for monomial, coefficient in poly.items():
        term = coefficient
        for value, power in zip(values, monomial):
            term *= value**power
        result += term
    return result


def literal_gradient_fixture(five):
    rows = ((1, 1), (2, 1), (1, 2), (3, 2), (2, 3))
    values = tuple(value for row in rows for value in row)
    first, second = five[0][1], five[1][1]
    q = add(scale(first, evaluate(second, values)),
            scale(second, -evaluate(first, values)))
    assert q and evaluate(q, values) == 0
    lam = []
    for i, (x, y) in enumerate(rows):
        qx = evaluate(derivative(q, 2 * i), values)
        qy = evaluate(derivative(q, 2 * i + 1), values)
        assert qx * x + qy * y == 0 and qy % x == 0
        value = qy // x
        assert qx == -value * y
        c, b, a = quadratic_coefficients(q, i)
        other_values = values[:2 * i] + values[2 * i + 2 :]
        disc = evaluate(b, other_values)**2 - 4 * evaluate(a, other_values) * evaluate(c, other_values)
        assert value**2 == disc
        lam.append(value)
    assert all(lam)
    assert sum(lam[i] * rows[i][0]**2 for i in range(5)) == 0
    assert sum(lam[i] * rows[i][0] * rows[i][1] for i in range(5)) == 0
    assert sum(lam[i] * rows[i][1]**2 for i in range(5)) == 0
    return tuple(lam)


def main():
    five = independent_graphs(5, 2)
    four = independent_graphs(4, 4)
    assert len(five) == 6 and len(four) == 5
    coefficients = [quadratic_coefficients(poly, 0) for _, poly in five]
    monomials = independent_monomials(five)
    coordinate_matrix = [[poly.get(m, 0) for _, poly in five]
                         for m in monomials]
    integer_index_bound = abs(determinant(coordinate_matrix))
    assert integer_index_bound.denominator == 1
    inverse_matrix = inverse(coordinate_matrix)
    conversion_constant = max(sum(abs(inverse_matrix[i][j]) for i in range(6))
                              for j in range(6))
    assert all(value.denominator == 1 for row in inverse_matrix for value in row)
    denominator = 1
    tests = 0
    for j, k in combinations_with_replacement(range(6), 2):
        polar = polarized_discriminant(coefficients[j], coefficients[k], j == k)
        expansion = solve_coefficients(four, polar)
        for value in expansion:
            denominator = lcm(denominator, value.denominator)
        tests += 1
    assert [internal_minimum(4, a) for a in range(5)] == [0, 0, 0, 4, 8]
    assert incident_minima() == {
        (False, 0): (0, 0), (False, 1): (0, 0),
        (False, 2): (0, 0), (False, 3): (2, 1),
        (False, 4): (4, 3), (True, 0): (-1, 0),
        (True, 1): (-1, 0), (True, 2): (0, 0),
        (True, 3): (2, 2), (True, 4): (4, 4),
    }
    assert literal_gradient_fixture(five) == (-4500, 500, -1625, -75, 1200)
    print(f"PASS: six five-row and five four-row graph basis elements; "
          f"five-row integral coordinate denominator {integer_index_bound}; "
          f"graph coefficient l1 conversion at most {conversion_constant} C_Q; "
          f"{tests} polarized discriminants, common denominator {denominator}; "
          "four-row internal minima 0,0,0,4,8; ten literal derivative minima; "
          "nonzero primitive-row gradient/discriminant/stress fixture.")


if __name__ == "__main__":
    main()
