"""Exact coherent boundary model in the actual labelled invariant system."""

from fractions import Fraction
from itertools import combinations
from math import comb


def add(a, b, scale=1):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += scale * x
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def multiply(a, b):
    if not a or not b:
        return ()
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def edge(i, j, inside, nodes):
    if i in inside and j in inside:
        return ()
    if i in inside:
        return (-1,)
    if j in inside:
        return (1,)
    return add(nodes[j], nodes[i], -1)


def graph(edges, inside, nodes):
    result = (1,)
    for i, j in edges:
        result = multiply(result, edge(i, j, inside, nodes))
    return result


def quartet(labels, inside, nodes):
    a, b, c, d = labels
    P = graph([(a, b), (c, d)], inside, nodes)
    R = graph([(a, c), (b, d)], inside, nodes)
    # At the actual directions z_i=i+1, both consecutive quartets used here
    # have P=1, R=4. This polynomial is 16 P^2-R^2.
    return add(tuple(16 * x for x in multiply(P, P)), multiply(R, R), -1)


def triangle_edges(labels):
    a, b, c = labels
    return [(a, b), (b, c), (c, a)]


def triangle_graph(labels, count):
    edges = []
    for k in range(count):
        edges.extend(triangle_edges(labels[3 * k:3 * k + 3]))
    for i in range(3 * count, len(labels), 2):
        edges.extend([(labels[i], labels[i + 1])] * 2)
    return edges


def image_coordinates(outside, nodes):
    last = outside[-1]
    base = multiply(add(nodes[outside[0]], nodes[last], -1),
                    add(nodes[outside[1]], nodes[last], -1))
    return [add(multiply(add(nodes[i], nodes[last], -1),
                         add(nodes[j], nodes[last], -1)), base, -1)
            for i, j in combinations(outside[:-1], 2)
            if (i, j) != (outside[0], outside[1])]


def rank(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    r = 0
    for j in range(len(a[0])):
        pivot = next((k for k in range(r, len(a)) if a[k][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        divisor = a[r][j]
        a[r] = [x / divisor for x in a[r]]
        for k in range(r + 1, len(a)):
            if a[k][j]:
                c = a[k][j]
                a[k] = [x - c * y for x, y in zip(a[k], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def matching_rank(r):
    pairs = list(combinations(range(r), 2))
    rows = []
    for a, b, c, d in combinations(range(r), 4):
        for e, f in [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]:
            row = [0] * len(pairs)
            for i, sign_i in [(e[0], -1), (e[1], 1)]:
                for j, sign_j in [(f[0], -1), (f[1], 1)]:
                    row[pairs.index(tuple(sorted((i, j))))] += sign_i * sign_j
            rows.append(row)
    return rank(rows)


def invariant_dimension(m):
    coefficients = [1]
    for _ in range(m):
        out = [0] * (len(coefficients) + 2)
        for i, x in enumerate(coefficients):
            out[i] += x
            out[i + 1] += x
            out[i + 2] += x
        coefficients = out
    return coefficients[m] - coefficients[m - 1]


def main():
    m, q = 16, 8
    A, A_prime = tuple(range(4)), tuple(range(4, 8))
    B = tuple(range(4, m))
    B_prime = tuple(i for i in range(m) if i not in A_prime)
    nodes = [(i + 1,) if i in A else (0, i + 1) for i in range(m)]
    actual = [(i + 1,) for i in range(m)]
    assert quartet(A, set(), actual) == quartet(A_prime, set(), actual) == ()
    first_graph = triangle_graph(B, 2)
    second_graph = triangle_graph(B, 4)
    disjoint_graph = triangle_graph(B_prime, 4)
    seen_orders, nonzero_restrictions = set(), 0
    for cut in combinations(range(m), q - 1):
        inside = set(cut)
        outside = tuple(i for i in range(m) if i not in inside)
        coordinates = image_coordinates(outside, nodes)
        order = min(i for value in coordinates for i, c in enumerate(value) if c)
        assert order in (0, 1, 2)
        seen_orders.add(order)
        Q = quartet(A, inside, nodes)
        relation = multiply(Q, graph(first_graph, inside, nodes))
        assert multiply(Q, graph(second_graph, inside, nodes)) == ()
        assert multiply(quartet(A_prime, inside, nodes),
                        graph(disjoint_graph, inside, nodes)) == ()
        a = len(inside.intersection(A))
        if a != 2:
            assert relation == ()
        else:
            assert order == 0
            if relation:
                nonzero_restrictions += 1
                assert len(relation) <= 3 and relation[0] == relation[1] == 0
        assert len(relation) <= order or relation[order] == 0

    assert seen_orders == {0, 1, 2} and nonzero_restrictions > 0
    # A first-smaller restriction of the complementary K1 space maps onto all
    # matching quadrics in seven variables; its exact rank is fourteen.
    assert matching_rank(7) == 14
    # Gcd-one graph certificates on both twelve-row complements: every
    # irreducible bracket is absent from at least one cyclic permutation.
    for complement in [B, B_prime]:
        edge_sets = []
        for shift in range(len(complement)):
            labels = complement[shift:] + complement[:shift]
            edge_sets.append({tuple(sorted(e)) for e in triangle_graph(labels, 4)})
        assert not set.intersection(*edge_sets)
    for m in [24, 30, 400]:
        q, n = m // 2, m // 2 + 1
        r = invariant_dimension(m)
        s = comb(m, q) // 2
        t = n * (n - 3) // 2
        lower = r - s - comb(m, q - 1) * t - 1
        assert 2 * lower > r
        if m == 24:
            assert (r, lower) == (834086421, 670484982)
        if m == 30:
            assert (r, lower) == (439742222071, 424540705110)
    # Interior-coherent existence model at ten rows, where K2=0 by the
    # hierarchy termination floor(m/6)+1=2. No matrix-rank or short-basis
    # assertion is inferred beyond the number of imposed linear equations.
    r10 = invariant_dimension(10)
    k1_10 = r10 - comb(10, 5) // 2
    cuts10 = comb(10, 4)
    quartet_kernel = invariant_dimension(4) - 1
    k1_6 = invariant_dimension(6) - comb(6, 3) // 2
    assert (r10, k1_10, cuts10) == (603, 477, 210)
    assert (quartet_kernel, k1_6) == (2, 5)
    assert k1_10 - cuts10 - 1 == 266 > quartet_kernel * k1_6 == 10
    interior_nodes = [(i + 1,) for i in range(10)]
    for cut in combinations(range(10), 4):
        outside = tuple(i for i in range(10) if i not in cut)
        coordinates = image_coordinates(outside, interior_nodes)
        assert len(coordinates) == 9 and any(coordinates)
    print(f"verified all 11440 cuts, {nonzero_restrictions} nonzero fixture restrictions, "
          "three limit orders, matching rank14, graph gcd1, dimensions through m400, "
          "and the m10 interior-coherent dimension model")


if __name__ == "__main__":
    main()
