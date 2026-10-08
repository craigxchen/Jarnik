"""Exact examples and restrictions for minimal_support_invariant_descent.md."""

from fractions import Fraction
from itertools import combinations

from check_smaller_cut_relation_rank import add, mul, restrict_graph, scale, variable
from check_complement_reflection_invariant_relations import conjugate, determinant, gaussian_gcd, norm


CYCLE = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0))
COMPLEMENT = ((0, 2), (2, 4), (4, 1), (1, 3), (3, 0))
MULTIPLIER_TRIPLES = ((0, 1, 2), (0, 1, 3), (0, 1, 4),
                      (0, 2, 3), (0, 2, 4), (0, 3, 4))


def rank(polynomials):
    monomials = sorted(set().union(*(poly.keys() for poly in polynomials)))
    basis = []
    pivots = []
    for poly in polynomials:
        row = [Fraction(poly.get(monomial, 0)) for monomial in monomials]
        for pivot, existing in zip(pivots, basis):
            factor = row[pivot]
            row = [a - factor * b for a, b in zip(row, existing)]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            row = [value / row[pivot] for value in row]
            pivots.append(pivot)
            basis.append(row)
    return len(basis)


def multiplier_edges(triple, offset=0):
    a, b, c = triple
    d, e = (j for j in range(5) if j not in triple)
    return tuple((i + offset, j + offset) for i, j in ((a, b), (b, c), (c, a), (d, e), (d, e)))


def graph_value(edges, rows):
    value = 1
    for i, j in edges:
        value *= determinant(rows[i], rows[j])
    return value


def g_restriction(inside, n=10):
    return add(scale(restrict_graph(CYCLE, inside, n), 12), restrict_graph(COMPLEMENT, inside, n))


def main():
    assert {frozenset(edge) for edge in CYCLE + COMPLEMENT} == set(map(frozenset, combinations(range(5), 2)))
    # Every possible bracket collision leaves exactly the other cycle nonzero.
    for i, j in combinations(range(5), 2):
        values = list(range(5))
        values[j] = values[i]
        points = [(value, 1) for value in values]
        first = graph_value(CYCLE, points)
        second = graph_value(COMPLEMENT, points)
        assert (first == 0) != (second == 0)
        assert 12 * first + second != 0

    one = {(0, 0): 1}
    a, b = variable(2, 0), variable(2, 1)
    chart_rows = [(one, {}), ({}, one), (one, one), (a, one), (b, one)]

    def chart_graph(edges):
        value = one
        for i, j in edges:
            bracket = add(mul(chart_rows[i][0], chart_rows[j][1]),
                          scale(mul(chart_rows[i][1], chart_rows[j][0]), -1))
            value = mul(value, bracket)
        return value

    chart_g = add(scale(chart_graph(CYCLE), 12), chart_graph(COMPLEMENT))
    assert chart_g == {(2, 0): -12, (1, 0): 12, (1, 1): 13, (1, 2): -1, (0, 1): -12}

    rows = [(1, 0), (0, 1)] + [(2*j, 1) for j in range(1, 9)]
    assert len(rows) == 10
    assert all(norm(gaussian_gcd(row, conjugate(row))) == 1 for row in rows)
    assert all(determinant(rows[i], rows[j]) for i, j in combinations(range(10), 2))
    assert 12 * graph_value(CYCLE, rows) + graph_value(COMPLEMENT, rows) == 0
    multipliers = [multiplier_edges(triple, 5) for triple in MULTIPLIER_TRIPLES]
    assert rank([restrict_graph(edges, set(), 10) for edges in multipliers]) == 6
    assert all(graph_value(edges, rows) for edges in multipliers)

    for inside in combinations(range(10), 5):
        g = g_restriction(set(inside))
        assert all(not mul(g, restrict_graph(edges, set(inside), 10)) for edges in multipliers)
    max_rank = 0
    for inside in combinations(range(10), 4):
        g = g_restriction(set(inside))
        restricted = [mul(g, restrict_graph(edges, set(inside), 10)) for edges in multipliers]
        local_rank = rank(restricted)
        assert local_rank <= 2
        max_rank = max(max_rank, local_rank)
    assert max_rank == 2

    # Projection groups exactly the prescribed number of original cut labels.
    for m in range(4, 10):
        for s in range(2, m + 1):
            groups = {mask: [] for mask in range(1 << s)}
            for old_mask in range(1 << m):
                groups[old_mask & ((1 << s) - 1)].append(old_mask)
            assert all(len(group) == 1 << (m - s) for group in groups.values())
            nonempty_old_groups = {new: [old for old in group if old] for new, group in groups.items()}
            assert all(len(group) == 1 << (m - s) for new, group in nonempty_old_groups.items() if new)
            assert len(nonempty_old_groups[0]) == (1 << (m - s)) - 1
            for i in range(s):
                regrouped = {old for new, group in groups.items() if new >> i & 1 for old in group}
                assert regrouped == {old for old in range(1 << m) if old >> i & 1}

    print("Five-cycle formula, all ten nondivisibility collisions, and actual primitive Gaussian zero pass.")
    print("Multiplier dimension6; all252 balanced restrictions vanish; all210 smaller images have rank≤2.")
    print("Exact support projection counts and row-factor preservation pass.")
    print("No asymptotic full-profile construction or uniform-bound counterexample is asserted.")


if __name__ == "__main__":
    main()
