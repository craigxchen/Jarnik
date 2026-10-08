"""Bounded exact search for low all-edge integer-cotangent objectives."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, lcm

from check_mobius_conductor_transfer import least_norm


def divisors(value):
    result = []
    for divisor in range(1, isqrt(value) + 1):
        if value % divisor:
            continue
        result.append(divisor)
        if divisor * divisor != value:
            result.append(value // divisor)
    return result


def quotient(x, y, scale):
    numerator = x * y + scale * scale
    denominator = x - y
    assert denominator and numerator % denominator == 0
    return numerator // denominator


def edge_norm(value, scale):
    common = gcd(abs(value), scale)
    a, b = value // common, scale // common
    divisor = 2 if a % 2 and b % 2 else 1
    assert (a * a + b * b) % divisor == 0
    return (a * a + b * b) // divisor


def all_edge_radius(xs, scale):
    values = [edge_norm(x, scale) for x in xs]
    values.extend(edge_norm(quotient(x, y, scale), scale)
                  for x, y in combinations(xs, 2))
    return lcm(*values)


def verify_clique(xs, scale, expected_radius=None):
    assert tuple(sorted(xs)) == tuple(xs)
    assert len(set(xs)) == len(xs)
    assert all(quotient(x, y, scale) for x, y in combinations(xs, 2))
    radius = all_edge_radius(xs, scale)
    if expected_radius is not None:
        assert radius == expected_radius
    return radius


def divisor_graph(scale, bound):
    vertices = list(range(1, bound + 1))
    adjacency = {vertex: set() for vertex in vertices}
    for lower in vertices:
        for gap in divisors(lower * lower + scale * scale):
            upper = lower + gap
            if upper <= bound:
                adjacency[lower].add(upper)
                adjacency[upper].add(lower)
    return vertices, adjacency


def search_sizes(scale=6, bound=10000, sizes=(4, 5)):
    vertices, adjacency = divisor_graph(scale, bound)
    best = {size: None for size in sizes}
    cliques = {size: 0 for size in sizes}
    maximum = max(sizes)

    def visit(current, candidates):
        size = len(current)
        if size in sizes:
            cliques[size] += 1
            radius = all_edge_radius(current, scale)
            objective = Fraction(radius * scale**4, current[0]**4)
            record = (objective, tuple(current), radius)
            if best[size] is None or record < best[size]:
                best[size] = record
        if size == maximum:
            return
        next_size = min((value for value in sizes if value > size), default=maximum)
        if len(candidates) < next_size - size:
            return
        for index, vertex in enumerate(candidates):
            next_candidates = [other for other in candidates[index + 1:]
                               if other in adjacency[vertex]]
            visit(current + [vertex], next_candidates)

    visit([], vertices)
    assert all(value is not None for value in best.values())
    return best, cliques


def fixed_fixture_audit():
    four = (1146, 1257, 1842, 4104)
    radius = verify_clique(four, 6, expected_radius=2163838625)
    assert radius == 5**3 * 13 * 17 * 29 * 37 * 73
    assert least_norm([(1, 0)] + [(x, 6) for x in four]) == radius
    objective = Fraction(radius * 6**4, four[0]**4)
    assert objective == Fraction(2163838625, 1330863361)
    assert [quotient(x, y, 6) for x, y in combinations(four, 2)] == [
        -12978, -3033, -1590, -3958, -1812, -3342]
    anchor = four[0]
    S = anchor * anchor + 6 * 6
    offsets = [x - anchor for x in four[1:]]
    cofactors = [S // gap for gap in offsets]
    assert offsets == [111, 696, 2958]
    assert cofactors == [11832, 1887, 444]
    assert all(S % gap == 0 for gap in offsets)
    incident = [2 * anchor + gap + cofactor
                for gap, cofactor in zip(offsets, cofactors)]
    normalized_gaps = [abs(x - y) // gcd(x, y)
                       for x, y in combinations(offsets, 2)]
    J = Fraction(incident[0] * incident[1] * incident[2],
                 normalized_gaps[0] * normalized_gaps[1] * normalized_gaps[2])
    assert J == 164250
    reanchored = (1146, 1590, 3033, 12978)
    assert tuple(sorted([anchor] + [anchor + value for value in cofactors])) == reanchored
    assert verify_clique(reanchored, 6, expected_radius=radius) == radius

    five = (9, 12, 18, 22, 48)
    radius5 = verify_clique(five, 6, expected_radius=325)
    assert least_norm([(1, 0)] + [(x, 6) for x in five]) == radius5
    assert Fraction(radius5 * 6**4, five[0]**4) == Fraction(5200, 81)
    return objective, radius5


if __name__ == "__main__":
    objective, radius5 = fixed_fixture_audit()
    best, cliques = search_sizes()
    assert best[4] == (Fraction(2163838625, 1330863361),
                       (1146, 1257, 1842, 4104), 2163838625)
    assert best[5] == (Fraction(5200, 81), (9, 12, 18, 22, 48), 325)
    print(f"PASS: all-edge fixtures, k=4 objective {objective}, k=5 radius {radius5}.")
    print(f"PASS: bounded L=6, X<=10000 search checked {cliques[4]} four-cliques "
          f"and {cliques[5]} five-cliques.")
    print("No global minimization or endpoint counterexample is asserted.")
