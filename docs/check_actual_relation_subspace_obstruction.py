"""Exact certificate for actual_relation_subspace_obstruction.md.

Standard library only.  This is an actual fixed rational configuration,
not an unbounded endpoint family or a full-profile existence assertion.
"""

from fractions import Fraction
from itertools import combinations
from math import comb, gcd

from check_smaller_cut_relation_rank import add, scale, restrict_graph, invariant_dimension


def rank(polynomials):
    monomials = sorted(set().union(*(set(p) for p in polynomials)))
    pivots = {}
    for poly in polynomials:
        row = [Fraction(poly.get(m, 0)) for m in monomials]
        for j, old in sorted(pivots.items()):
            if row[j]:
                multiple = row[j]
                row = [a - multiple * b for a, b in zip(row, old)]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            divisor = row[pivot]
            pivots[pivot] = [value / divisor for value in row]
    return len(pivots)


def evaluate(poly, points):
    return sum(coefficient * product(x ** power for x, power in zip(points, monomial))
               for monomial, coefficient in poly.items())


def product(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def primitive(poly):
    content = 0
    for coefficient in poly.values():
        content = gcd(content, abs(coefficient))
    assert content
    return {monomial: coefficient // content for monomial, coefficient in poly.items()}


def check_eight_rows():
    first_triangle = [(0, 1), (1, 2), (2, 0)]
    residual = list(range(3, 8))
    basis_graphs, basis_polynomials, labels = [], [], []
    for triple in combinations(residual, 3):
        a, b, c = triple
        d, e = [i for i in residual if i not in triple]
        graph = [(a, b), (b, c), (c, a), (d, e), (d, e)]
        poly = restrict_graph(graph, set(), n=8)
        if rank(basis_polynomials + [poly]) > len(basis_polynomials):
            basis_graphs.append(graph)
            basis_polynomials.append(poly)
            labels.append(tuple(i + 1 for i in triple))
    assert len(basis_graphs) == invariant_dimension(5) == 6

    # P_i=4(i+1)+i are conjugate-primitive, odd norm, distinct directions.
    points = [4 * (i + 1) for i in range(8)]
    values = [evaluate(poly, points) for poly in basis_polynomials]
    assert all(values)
    full_graphs = [first_triangle + graph for graph in basis_graphs]
    full_polynomials = [restrict_graph(graph, set(), n=8) for graph in full_graphs]
    raw_relations = [add(scale(full_polynomials[j], values[0]),
                         scale(full_polynomials[0], -values[j])) for j in range(1, 6)]
    relations = [primitive(poly) for poly in raw_relations]
    assert rank(relations) == 5
    assert all(evaluate(poly, points) == 0 for poly in relations)

    # Use unnormalized relations for the restriction ranks; scaling is harmless.
    balanced_count = 0
    for cut in combinations(range(8), 4):
        restrictions = [restrict_graph(graph, set(cut), n=8) for graph in full_graphs]
        assert all(not poly for poly in restrictions)
        balanced_count += 1
    ranks = {}
    for cut in combinations(range(8), 3):
        restrictions = [restrict_graph(graph, set(cut), n=8) for graph in full_graphs]
        images = [add(scale(restrictions[j], values[0]),
                      scale(restrictions[0], -values[j])) for j in range(1, 6)]
        actual_rank = rank(images)
        assert actual_rank <= 2
        ranks[actual_rank] = ranks.get(actual_rank, 0) + 1
    assert ranks == {0: 26, 2: 30}
    print("Residual triangle labels:", labels)
    print("Five primitive relation coefficient sums:", [sum(map(abs, p.values())) for p in relations])
    print("Actual zero relations: rank 5; balanced cuts:", balanced_count)
    print("All 56 first-smaller restriction ranks:", ranks)


def height_coefficient(n):
    return n * 2 ** (n - 2) - sum(comb(n, a) * max(0, 2 * a - n)
                                 for a in range(n + 1))


def check_thresholds():
    for n in range(3, 61, 2):
        loss = sum(comb(n, a) * max(0, 2 * a - n) for a in range(n + 1))
        assert loss == n * comb(n - 1, (n - 1) // 2)
    first = None
    print("m, n, residual dimension, A_n, t, 8*t*A_n/(r_n-t):")
    for m in range(8, 61, 2):
        n, q = m - 3, m // 2
        r, height = invariant_dimension(n), height_coefficient(n)
        t = (q + 1) * (q - 2) // 2
        assert r > t
        ratio = Fraction(8 * t * height, r - t)
        if ratio < 2 and first is None:
            first = m
        if m in (8, 24, 30, 36, 38, 40, 42, 44, 48):
            print(m, n, r, height, t, f"{float(ratio):.12f}")
    assert first == 40
    assert 8 * 189 * 935530313516 < 2 * (707854577312178 - 189)
    print("First even m in this exact comparison with strict ratio below 2:", first)


if __name__ == "__main__":
    check_eight_rows()
    check_thresholds()
    print("No unbounded endpoint family or existence of a large full profile is asserted.")
