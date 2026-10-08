#!/usr/bin/env python3
"""Verify two exact K5 zero-edge/opposite-disjoint-edge identities.

Uses only integer polynomial arithmetic from the Python standard library.
No symbolic solver or numerical feasibility result is trusted.
"""
import itertools
import json
from pathlib import Path

ZERO = (0, 0, 0)
EDGES = list(itertools.combinations(range(5), 2))


def add(*polys):
    result = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            result[monomial] = result.get(monomial, 0) + coefficient
    return {m: c for m, c in result.items() if c}


def scale(poly, coefficient):
    return {m: c * coefficient for m, c in poly.items() if c * coefficient}


def mul(left, right):
    result = {}
    for a, c in left.items():
        for b, d in right.items():
            m = tuple(x + y for x, y in zip(a, b))
            result[m] = result.get(m, 0) + c * d
    return {m: c for m, c in result.items() if c}


def decode(terms):
    result = {}
    for a, b, c, coefficient in terms:
        assert all(isinstance(v, int) for v in (a, b, c, coefficient))
        assert min(a, b, c) >= 0 and (a, b, c) not in result
        result[a, b, c] = coefficient
    return {m: c for m, c in result.items() if c}


def main():
    cases = json.loads(Path(__file__).with_name(
        'repeated_k5_zero_edge_certificates.json').read_text())
    assert len(cases) == 2
    for index, case in enumerate(cases):
        edges = [decode(p) for p in case['edge_forms']]
        multipliers = [decode(p) for p in case['multipliers']]
        assert len(edges) == 10 and len(multipliers) == 4
        expected_pairs = [((0, 2), (1, 3)), ((1, 2), (3, 4))]
        pair = tuple(map(tuple, case['opposite_edges']))
        assert pair == expected_pairs[index]
        assert edges[0] == {}
        assert edges[EDGES.index(pair[0])] == {ZERO: 1}
        assert edges[EDGES.index(pair[1])] == {ZERO: -1}
        # Literal free coordinates establish that the three-parameter
        # affine solution is injective. Rank seven below proves completeness.
        free_indices = (7, 8, 9) if index == 0 else (6, 7, 8)
        for j, edge_index in enumerate(free_indices):
            assert edges[edge_index] == {
                tuple(int(k == j) for k in range(3)): 1}
        vertex_sums = [add(*(edges[k] for k, e in enumerate(EDGES)
                              if i in e)) for i in range(5)]
        assert all(p == vertex_sums[0] for p in vertex_sums)
        cubes = [mul(mul(p, p), p) for p in edges]
        vertex_cubes = [add(*(cubes[k] for k, e in enumerate(EDGES)
                               if i in e)) for i in range(5)]
        differences = [add(p, scale(vertex_cubes[0], -1))
                       for p in vertex_cubes[1:]]
        assert all(c % 3 == 0 for p in differences for c in p.values())
        reduced = [{m: c // 3 for m, c in p.items()} for p in differences]
        actual = add(*(mul(a, b) for a, b in zip(multipliers, reduced)))
        target = decode(case['target'])
        expected_target = {(1, 0, 0): 3200} if index == 0 else {ZERO: 9600}
        assert target == expected_target and actual == target
        # Exact rank of the seven independent affine constraints.
        from fractions import Fraction
        matrix = [[Fraction(int(i in e) - int(0 in e)) for e in EDGES]
                  for i in range(1, 5)]
        matrix += [[Fraction(k == j) for k in range(10)]
                   for j in (0, EDGES.index(pair[0]), EDGES.index(pair[1]))]
        rank = 0
        for col in range(10):
            pivot = next((r for r in range(rank, len(matrix))
                          if matrix[r][col]), None)
            if pivot is None:
                continue
            matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
            divisor = matrix[rank][col]
            matrix[rank] = [v / divisor for v in matrix[rank]]
            for r in range(rank + 1, len(matrix)):
                factor = matrix[r][col]
                matrix[r] = [a - factor * b
                             for a, b in zip(matrix[r], matrix[rank])]
            rank += 1
        assert rank == 7
        print(f'Case {index + 1}: exact degree-seven identity verified; '
              f'target = {expected_target}.')
    # Every disjoint pair avoiding edge 01 falls into exactly these two
    # types: its four endpoints either include both 0,1 or exactly one.
    pairs = [(a, b) for a, b in itertools.combinations(EDGES[1:], 2)
             if set(a).isdisjoint(b)]
    assert {len(set((0, 1)) & (set(a) | set(b)))
            for a, b in pairs} == {1, 2}
    permutations = [p for p in itertools.permutations(range(5))
                    if set(p[:2]) == {0, 1}]
    orbit_sets = []
    for representative in (((0, 2), (1, 3)), ((1, 2), (3, 4))):
        orbit = {tuple(sorted(tuple(sorted((p[a], p[b])))
                              for a, b in representative)) for p in permutations}
        assert len(orbit) == 6
        orbit_sets.append(orbit)
    assert orbit_sets[0].isdisjoint(orbit_sets[1])
    assert orbit_sets[0] | orbit_sets[1] == set(pairs)
    print('Both six-element zero-edge/disjoint-pair orbits are covered.')


if __name__ == '__main__':
    main()
