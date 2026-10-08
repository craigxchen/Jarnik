"""Exact algebra accompanying odd_support_relation_rigidity.md; stdlib only."""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import comb

from check_smaller_cut_relation_rank import restrict_graph


def rank(rows):
    pivots = {}
    for original in rows:
        row = [Fraction(value) for value in original]
        for pivot, previous in sorted(pivots.items()):
            scalar = row[pivot]
            if scalar:
                row = [a - scalar * b for a, b in zip(row, previous)]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            scalar = row[pivot]
            pivots[pivot] = [value / scalar for value in row]
    return len(pivots)


def polynomial_rank(polynomials):
    monomials = sorted(set().union(*(poly.keys() for poly in polynomials)))
    return rank([[poly.get(monomial, 0) for monomial in monomials]
                 for poly in polynomials])


def graphs(degrees):
    """All loopless multigraphs with the specified labelled degrees."""
    if not any(degrees):
        yield ()
        return
    i = next(j for j, degree in enumerate(degrees) if degree)
    available = [j for j in range(i + 1, len(degrees)) if degrees[j]]
    for targets in combinations_with_replacement(available, degrees[i]):
        remaining = list(degrees)
        remaining[i] = 0
        for j in targets:
            remaining[j] -= 1
        if min(remaining) >= 0:
            for rest in graphs(tuple(remaining)):
                yield tuple((i, j) for j in targets) + rest


def invariant_dimension(degrees):
    coefficients = [1]
    for degree in degrees:
        following = [0] * (len(coefficients) + degree)
        for i, coefficient in enumerate(coefficients):
            for j in range(degree + 1):
                following[i + j] += coefficient
        coefficients = following
    weight = sum(degrees) // 2
    return coefficients[weight] - coefficients[weight - 1]


def balanced_value(edges, inside):
    value = 1
    for i, j in edges:
        if (i in inside) == (j in inside):
            return 0
        if j in inside:
            value = -value
    return value


def mixed_degree_checks():
    for degrees, dimension, scalar_rank in [
        ((2, 1, 1, 1, 1), 3, 3),
        ((2, 2, 2, 1, 1), 4, 3),
    ]:
        monomials = list(graphs(degrees))
        polys = [restrict_graph(edges, set(), 5) for edges in monomials]
        cuts = [set(j for j, bit in enumerate(bits) if bit)
                for bits in product((0, 1), repeat=5)
                if sum(bit * degree for bit, degree in zip(bits, degrees))
                == sum(degrees) // 2]
        matrix = [[balanced_value(edges, cut) for cut in cuts]
                  for edges in monomials]
        assert invariant_dimension(degrees) == dimension
        assert polynomial_rank(polys) == dimension
        assert rank(matrix) == scalar_rank
    kernel_graph = ((0, 1), (1, 2), (2, 0), (3, 4))
    assert restrict_graph(kernel_graph, set(), 5)
    for bits in product((0, 1), repeat=5):
        if sum(bit * degree for bit, degree in zip(bits, (2, 2, 2, 1, 1))) == 4:
            cut = {j for j, bit in enumerate(bits) if bit}
            assert balanced_value(kernel_graph, cut) == 0


def odd_hierarchy_checks():
    for m in (3, 5, 7):
        q = m // 2
        monomials = list(graphs((2,) * m))
        cuts = list(combinations(range(m), q))
        # At the first odd level all restrictions are homogeneous linear.
        matrix = []
        for edges in monomials:
            row = []
            for cut_tuple in cuts:
                cut = set(cut_tuple)
                poly = restrict_graph(edges, cut, m)
                assert all(sum(monomial) == 1 for monomial in poly)
                coefficients = [poly.get(tuple(int(i == j) for i in range(m)), 0)
                                for j in range(m)]
                assert sum(coefficients) == 0
                row.extend(coefficients)
            matrix.append(row)
        assert rank(matrix) == invariant_dimension((2,) * m)
    for m in range(3, 30, 2):
        q = m // 2
        ell = (m - 3) // 6
        r = (m - 3 - 6 * ell) // 2
        assert r in (0, 1, 2)
        edges = []
        selected = set()
        index = 0
        for _ in range(2 * ell + 1):
            edges.extend(((index, index + 1), (index + 1, index + 2),
                          (index + 2, index)))
            selected.add(index)
            index += 3
        for _ in range(r):
            edges.extend(((index, index + 1),) * 2)
            selected.add(index)
            index += 2
        assert index == m and len(selected) == q - ell
        assert restrict_graph(edges, selected, m)
        for k in range(q + 1):
            n, d = q + k + 1, 2 * k + 1
            assert (2 * d <= n) == (3 * k + 1 <= q)
            assert 4 * n - 2 * d == 2 * m
            assert 2 * m + (2 * m) * (2 * m - 1) + 2 * m == 2 * m * (2 * m + 1)
            if 2 * d <= n:
                assert comb(n, d) - comb(n, d - 1) > 0


def five_row_geometry_and_capacity():
    boundary = list(combinations(range(5), 2))
    adjacency = {(a, b) for a, b in combinations(range(10), 2)
                 if set(boundary[a]).isdisjoint(boundary[b])}
    assert len(adjacency) == 15
    assert all(not all(tuple(sorted(pair)) in adjacency
                       for pair in combinations(triple, 2))
               for triple in combinations(range(10), 3))
    # On Bl_3(P1 x P1), -K=(2,2; -1,-1,-1) has square5.
    assert 2 * 2 * 2 - 3 == 5
    assert len(boundary) // 2 + 1 > 5
    assert 2 * 5 * (2 * 5 + 1) == 110
    assert 3**5 - 1 == 242

    cycle = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0))
    triangle = ((5, 6), (6, 7), (7, 5))
    for edges, m, capacity in [(cycle, 5, 2), (triangle, 8, 6),
                               (cycle + triangle, 8, 3)]:
        # The middle graph leaves five zero-degree labels unused, adding5.
        maximum = max(len(cut) for length in range(m + 1)
                      for cut in combinations(range(m), length)
                      if restrict_graph(edges, set(cut), m))
        assert maximum == capacity
    # Actual distinct rational directions, not an asserted full profile:
    # the prior fixed irreducible relation has a nonzero two-cut restriction.
    from check_minimal_support_invariant_descent import g_restriction
    from check_smaller_cut_relation_rank import mul
    nonzero_cut = next(set(cut) for cut in combinations(range(5), 2)
                       if g_restriction(set(cut), 8))
    full_cut = nonzero_cut | {5}
    assert mul(g_restriction(full_cut, 8), restrict_graph(triangle, full_cut, 8))
    assert not any(mul(g_restriction(set(cut), 8), restrict_graph(triangle, set(cut), 8))
                   for cut in combinations(range(8), 4))


def main():
    mixed_degree_checks()
    print('Five-row mixed-degree invariant/scalar ranks (3,3) and (4,3), and the nonzero kernel monomial pass.')
    odd_hierarchy_checks()
    print('Odd first-level injectivity for 3, 5, 7 rows and sharp termination witnesses through 29 rows pass.')
    five_row_geometry_and_capacity()
    print('Boundary incidence budget, finite constants, and exact triangle multiplication capacity pass.')
    print('The geometric pair theorem is proved in the note; no individual small relation or uniform arc bound is asserted.')


if __name__ == '__main__':
    main()
