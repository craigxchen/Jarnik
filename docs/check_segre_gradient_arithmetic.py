"""Exact polynomial certificate for the six-row Segre/gradient coordinates."""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import comb, lcm

from check_near_balanced_invariant_congruences import (
    PARTITIONS, ZERO, add, mul, scale, variable,
)


MATCHINGS = (
    ((0, 1), (2, 3), (4, 5)),
    ((0, 1), (2, 5), (3, 4)),
    ((0, 3), (1, 2), (4, 5)),
    ((0, 5), (1, 2), (3, 4)),
    ((0, 5), (1, 4), (2, 3)),
)


def product(factors):
    result = {ZERO: 1}
    for factor in factors:
        result = mul(result, factor)
    return result


def edge(i, j):
    return add(variable(i), scale(variable(j), -1))


def graph(edges):
    return product(edge(i, j) for i, j in edges)


def rref(rows, coefficient_columns=None):
    matrix = [[Fraction(value) for value in row] for row in rows]
    ncols = len(matrix[0]) if matrix else 0
    if coefficient_columns is None:
        coefficient_columns = ncols
    pivot_row = 0
    pivots = []
    for column in range(coefficient_columns):
        selected = next((i for i in range(pivot_row, len(matrix)) if matrix[i][column]), None)
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = matrix[selected], matrix[pivot_row]
        divisor = matrix[pivot_row][column]
        matrix[pivot_row] = [value / divisor for value in matrix[pivot_row]]
        for i in range(len(matrix)):
            if i != pivot_row and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [x - factor * y for x, y in zip(matrix[i], matrix[pivot_row])]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return matrix, pivots


def coordinates(basis, targets):
    monomials = sorted(set().union(*(set(p) for p in basis + targets)))
    rows = [[p.get(m, 0) for p in basis + targets] for m in monomials]
    reduced, pivots = rref(rows, len(basis))
    assert len(pivots) == len(basis)
    assert not any(any(row[len(basis):]) for row in reduced[len(basis):])
    answer = [tuple(reduced[i][len(basis) + j] for i in range(len(basis)))
              for j in range(len(targets))]
    for target, coefficients in zip(targets, answer):
        reconstructed = {}
        for coefficient, polynomial in zip(coefficients, basis):
            reconstructed = add(reconstructed, scale(polynomial, coefficient))
        assert reconstructed == target
    return answer


def rank(polynomials):
    monomials = sorted(set().union(*(set(p) for p in polynomials)))
    rows = [[p.get(m, 0) for p in polynomials] for m in monomials]
    return len(rref(rows)[1])


def all_matchings(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for j in vertices[1:]:
        rest = tuple(i for i in vertices if i not in (first, j))
        for matching in all_matchings(rest):
            yield ((first, j),) + matching


def matching_value(matching, inside):
    vectors = [(1, 0) if i in inside else (0, 1) for i in range(6)]
    result = 1
    for i, j in matching:
        result *= vectors[i][0] * vectors[j][1] - vectors[j][0] * vectors[i][1]
    return result


def main():
    matching_polynomials = [graph(matching) for matching in MATCHINGS]
    A, B, C, D, E = matching_polynomials
    total = {}
    for polynomial in matching_polynomials:
        total = add(total, polynomial)
    cubic = add(product((B, C, E)), scale(product((A, D, total)), -1))
    assert not cubic

    gradients = [
        scale(mul(D, add(total, A)), -1),
        add(mul(C, E), scale(mul(A, D), -1)),
        add(mul(B, E), scale(mul(A, D), -1)),
        scale(mul(A, add(total, D)), -1),
        add(mul(B, C), scale(mul(A, D), -1)),
    ]
    p, u, v, q, w = gradients
    z = add(mul(p, q), scale(add(add(mul(u, v), mul(u, w)), mul(v, w)), -1))
    t = product((u, v, w))
    h = add(add(add(u, v), w), scale(add(p, q), -1))
    inverse_gradients = [
        add(scale(mul(z, q), 2), scale(t, 4)),
        add(add(scale(mul(z, add(v, w)), -2), scale(product((v, w, h)), -4)), scale(t, -4)),
        add(add(scale(mul(z, add(u, w)), -2), scale(product((u, w, h)), -4)), scale(t, -4)),
        add(scale(mul(z, p), 2), scale(t, 4)),
        add(add(scale(mul(z, add(u, v)), -2), scale(product((u, v, h)), -4)), scale(t, -4)),
    ]
    vandermonde = graph(tuple(combinations(range(6), 2)))
    for inverse, matching in zip(inverse_gradients, matching_polynomials):
        assert inverse == scale(mul(vandermonde, matching), -4)
    # Euler's identities now also certify J(grad Phi)=0, without expanding
    # the larger quartic composition: grad J . grad Phi = -4 V M . grad Phi.
    triangles = []
    for first, second in PARTITIONS:
        edges = []
        for a, b, c in (first, second):
            edges.extend(((a, b), (b, c), (c, a)))
        triangles.append(graph(edges))
    triangle_coordinates = coordinates(gradients, triangles)
    # Select an exact triangle basis.
    triangle_basis = []
    basis_indices = []
    for index, polynomial in enumerate(triangles):
        if rank(triangle_basis + [polynomial]) > len(triangle_basis):
            triangle_basis.append(polynomial)
            basis_indices.append(index)
    assert len(triangle_basis) == 5
    gradient_coordinates = coordinates(triangle_basis, gradients)

    every_matching = list(all_matchings(tuple(range(6))))
    matching_coordinates = coordinates(matching_polynomials, [graph(m) for m in every_matching])
    assert all(c.denominator == 1 for row in matching_coordinates for c in row)

    quadratic_pairs = list(combinations_with_replacement(range(5), 2))
    quadratics = [mul(matching_polynomials[i], matching_polynomials[j]) for i, j in quadratic_pairs]
    assert rank(quadratics) == 15
    balanced_matrix = []
    for inside in combinations(range(6), 3):
        values = [matching_value(m, set(inside)) for m in MATCHINGS]
        balanced_matrix.append([values[i] * values[j] for i, j in quadratic_pairs])
        a, b, c, d, e = values
        assert (-d * (2*a+b+c+d+e), c*e-a*d, b*e-a*d,
                -a * (a+b+c+2*d+e), b*c-a*d) == (0,) * 5
    balanced_reduced, balanced_pivots = rref(balanced_matrix)
    assert len(balanced_pivots) == 10
    assert all(entry.denominator == 1 for row in balanced_reduced for entry in row)
    free_columns = [j for j in range(15) if j not in balanced_pivots]
    gradient_quadratic_coordinates = coordinates(quadratics, gradients)
    free_matrix = [[row[j] for j in free_columns] for row in gradient_quadratic_coordinates]
    assert all(sum(abs(value) for value in row) == 1 for row in free_matrix)
    assert all(sum(abs(free_matrix[i][j]) for i in range(5)) == 1 for j in range(5))

    f = (0, 0, 0, 1, 2, 4, 6)
    for size in range(7):
        for inside_tuple in combinations(range(6), size):
            inside = set(inside_tuple)
            matching_orders = [sum(i in inside and j in inside for i, j in matching)
                               for matching in every_matching]
            triangle_orders = [sum(comb(len(set(part) & inside), 2) for part in pair)
                               for pair in PARTITIONS]
            assert min(matching_orders) == max(0, size - 3)
            assert min(triangle_orders) == f[size]
    assert sum(comb(6, size) * max(0, size - 3) for size in range(7)) == 30
    assert sum(comb(6, size) * f[size] for size in range(7)) == 80
    assert [f[size] - 2 * max(0, size-3) for size in range(7)] == [0,0,0,1,0,0,0]

    print("Cubic identity passed: BCE - AD(A+B+C+D+E)=0.")
    print("Dual quartic inverse passed: grad J(grad Phi)=-4 Vandermonde * matching vector.")
    print("Five matching basis elements span all15 matchings integrally.")
    print("Quadratic rank15; balanced evaluation rank10; gradient/kernel rank5.")
    print("Triangle products in gradient coordinates, triple order 123,124,...,156:")
    for row in triangle_coordinates:
        print(tuple(str(c) for c in row))
    print("Selected triangle basis indices:", basis_indices)
    print("Gradients in that triangle basis:")
    for row in gradient_coordinates:
        print(tuple(str(c) for c in row))
    denominator = lcm(*(c.denominator for row in triangle_coordinates + gradient_coordinates for c in row))
    print("Common denominator for both span directions:", denominator)
    print("All64 cut minima passed; matching sum30, triangle sum80, gradient excess20.")


if __name__ == "__main__":
    main()
