"""Exact combinatorial and elliptic checks for five_row_del_pezzo_arithmetic.md."""

from fractions import Fraction
from itertools import combinations, permutations
from math import comb

from check_smaller_cut_relation_rank import add, mul, restrict_graph, scale, variable
from check_minimal_support_invariant_descent import rank
from check_two_small_invariant_relations import determinant


EDGES = tuple(combinations(range(5), 2))


def graphs():
    result = []
    for a, b, c in combinations(range(5), 3):
        d, e = (j for j in range(5) if j not in (a, b, c))
        result.append(((a, b), (b, c), (c, a), (d, e), (d, e)))
    for tail in permutations(range(1, 5)):
        if tail[0] < tail[-1]:
            cycle = (0,) + tail
            result.append(tuple((cycle[j], cycle[(j + 1) % 5]) for j in range(5)))
    return result


def internal(graph, cut):
    return sum(i in cut and j in cut for i, j in graph)


def power(poly, exponent):
    dimension = len(next(iter(poly), (0, 0, 0)))
    result = {(0,) * dimension: 1}
    for _ in range(exponent):
        result = mul(result, poly)
    return result


def elliptic_equation(x, y, z, t):
    one = {(0, 0, 0): 1}
    return add(add(add(mul(power(y, 2), z), mul(add(one, scale(t, -1)), mul(mul(x, y), z))),
                       scale(mul(mul(t, y), power(z, 2)), -1)),
               add(scale(power(x, 3), -1), mul(t, mul(power(x, 2), z))))


def main():
    graph_list = graphs()
    assert len(graph_list) == 22
    polynomials = [restrict_graph(graph, set(), 5) for graph in graph_list]
    assert rank(polynomials) == 6
    for edge in EDGES:
        selected = [poly for graph, poly in zip(graph_list, polynomials)
                    if frozenset(edge) in map(frozenset, graph)]
        assert rank(selected) == 4

    g0_sum = 0
    boundary_sums = {edge: 0 for edge in EDGES}
    first_cycle = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0))
    second_cycle = ((0, 2), (2, 4), (4, 1), (1, 3), (3, 0))
    for mask in range(32):
        cut = {i for i in range(5) if mask >> i & 1}
        g0 = max(0, 2 * len(cut) - 5)
        assert min(internal(graph, cut) for graph in graph_list) == g0
        assert min(internal(first_cycle, cut), internal(second_cycle, cut)) == g0
        g0_sum += g0
        for edge in EDGES:
            selected_graphs = [graph for graph in graph_list if frozenset(edge) in map(frozenset, graph)]
            expected = g0 + int(cut == set(edge) or cut == set(range(5)) - set(edge))
            assert min(internal(graph, cut) for graph in selected_graphs) == expected
            boundary_sums[edge] += expected
    assert g0_sum == 30 and set(boundary_sums.values()) == {32}

    cross_intersections = [(edge, other) for edge in first_cycle for other in second_cycle
                           if set(edge).isdisjoint(other)]
    assert len(cross_intersections) == 5
    assert all(sum(edge == first for first, _ in cross_intersections) == 1 for edge in first_cycle)
    assert all(sum(edge == second for _, second in cross_intersections) == 1 for edge in second_cycle)

    for d in range(13):
        choose = lambda a, b: comb(a, b) if a >= b else 0
        hilbert = choose(5*d+3, 3) - 5*choose(3*d+2, 3) + 10*choose(d+1, 3)
        assert hilbert == 5*d*(d+1)//2 + 1

    # Exact polynomial Tate transformation, with affine plane coordinates a,b.
    a, b, t = (variable(3, j) for j in range(3))
    one = {(0, 0, 0): 1}
    coefficient = add(add(t, mul(add(t, one), b)), scale(power(b, 2), -1))
    g = add(add(scale(mul(t, power(a, 2)), -1), mul(coefficient, a)), scale(mul(t, b), -1))
    denominator = add(b, scale(mul(t, add(a, scale(one, -1))), -1))
    x, y = mul(t, b), power(t, 2)
    assert elliptic_equation(x, y, denominator, t) == mul(power(t, 4), g)

    b2 = add(add(power(t, 2), scale(t, -6)), one)
    b4 = add(power(t, 2), scale(t, -1))
    b6, b8 = power(t, 2), scale(power(t, 3), -1)
    disc = add(add(scale(mul(power(b2, 2), b8), -1), scale(power(b4, 3), -8)),
               add(scale(power(b6, 2), -27), scale(mul(mul(b2, b4), b6), 9)))
    assert disc == mul(power(t, 5), add(add(power(t, 2), scale(t, -11)), scale(one, -1)))
    assert 12**5 * (12**2 - 11*12 - 1) != 0

    # Independently certify the four distinct branch points at t=12.
    quartic = [1, -26, 145, -264, 144]
    derivative = [4, -78, 290, -264]
    sylvester = [[0]*j + quartic + [0]*(2-j) for j in range(3)]
    sylvester += [[0]*j + derivative + [0]*(3-j) for j in range(4)]
    assert determinant(sylvester) == 700710912

    # The five rational boundary points on the Tate cubic, symbolically.
    zero = {}
    for px, py, pz in ((zero, one, zero), (zero, zero, one),
                       (t, power(t, 2), one), (t, zero, one), (zero, t, one)):
        assert not elliptic_equation(px, py, pz, t)
    point_images = (
        ((1, 1, 1), (t, power(t, 2), one)),
        ((1, 0, 0), (zero, zero, one)),
        ((1, 0, 1), (zero, one, zero)),
        ((0, 0, 1), (zero, t, one)),
        ((0, 1, 0), (t, zero, one)),
    )
    for (pa, pb, pz), target in point_images:
        transformed = (scale(t, pb), scale(power(t, 2), pz),
                       add(scale(one, pb), scale(t, pz-pa)))
        for i, j in combinations(range(3), 2):
            assert mul(transformed[i], target[j]) == mul(transformed[j], target[i])
    # The interior point from (a,b)=(2,3) is two-torsion at t=12.
    px, py = Fraction(-4), Fraction(-16)
    assert py*py - 11*px*py - 12*py == px**3 - 12*px*px
    assert 2*py - 11*px - 12 == 0

    print("22 graph coordinates span6; every boundary hyperplane space spans4.")
    print("All32 core minima and all320 boundary minima pass; sums30 and32.")
    print("Hilbert formula, five paired pentagon basepoints, Tate transform and discriminant pass.")
    print("G_12 is smooth genus one; its five boundary points are the stated5-torsion orbit.")
    print("No full-profile point family is asserted.")


if __name__ == "__main__":
    main()
