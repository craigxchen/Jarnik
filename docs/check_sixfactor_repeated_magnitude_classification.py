#!/usr/bin/env python3
"""Exhaustive rank-four equal-sum classification in {0,1,2} x {0,1}^4.

Only standard-library integer arithmetic; coefficients are not bounded.
"""
import itertools
import math
from collections import Counter, defaultdict

POINTS = list(itertools.product(range(3), range(2), range(2), range(2), range(2)))
NEIGHBORS = [{j for j in range(i + 1, len(POINTS))
              if sum(abs(a - b) for a, b in zip(POINTS[i], POINTS[j])) >= 3}
             for i in range(len(POINTS))]


def determinant(matrix):
    """Bareiss elimination with exact divisions and row pivoting."""
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a) - 1):
        if not a[k][k]:
            pivot = next((r for r in range(k + 1, len(a)) if a[r][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                numerator = a[i][j] * pivot_value - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        previous = pivot_value
    return sign * a[-1][-1]


def main():
    counts = Counter()
    number_cliques = 0

    def visit(indices, candidates):
        nonlocal number_cliques
        if len(indices) == 5:
            number_cliques += 1
            reference = POINTS[indices[0]]
            matrix = [[a - b for a, b in zip(POINTS[i], reference)]
                      for i in indices[1:]]
            normal = [(-1) ** k * determinant([row[:k] + row[k + 1:]
                                               for row in matrix])
                      for k in range(5)]
            assert all(sum(a * b for a, b in zip(row, normal)) == 0
                       for row in matrix)
            # Nonzero cofactor vector is equivalent to affine rank four.
            # A strictly positive normal and distinct coefficients are required.
            if not (all(a > 0 for a in normal) or all(a < 0 for a in normal)):
                return
            if len(set(normal)) != 5:
                return
            divisor = math.gcd(*normal)
            normal = tuple(abs(a) // divisor for a in normal)
            canonical = (normal[0], *sorted(normal[1:]))
            counts[canonical] += 1
            return
        for i in sorted(candidates):
            visit(indices + [i], candidates & NEIGHBORS[i])

    visit([], set(range(len(POINTS))))
    assert number_cliques == 3360
    assert counts == {(1, 2, 3, 4, 5): 48, (2, 1, 3, 4, 5): 48}
    print(f'{number_cliques} distance-three five-cliques; '
          f'{sum(counts.values())} admissible rank-four cliques.')
    subsets = [{4, 6}, {1, 3, 6}, {1, 4, 5}, {2, 3, 5}, {1, 2, 3, 4}]
    for alpha, (u, v), expected_sums in [
            ((1, 2, 3, 4, 5), (2, -5), {7, 9}),
            ((2, 1, 3, 4, 5), (1, 2), {8, 9})]:
        fibers = defaultdict(set)
        for point in POINTS:
            fibers[sum(a * b for a, b in zip(alpha, point))].add(point)
        large = {key: rows for key, rows in fibers.items() if len(rows) >= 5}
        assert set(large) == expected_sums
        assert all(len(rows) == 5 for rows in large.values())
        coefficients = (u, 2 * u, v, u + v, 2 * u + v, 3 * u + v)
        rows = set()
        for subset in subsets:
            bits = [int(j + 1 in subset) for j in range(6)]
            positive_bits = [bit if a > 0 else 1 - bit
                             for a, bit in zip(coefficients, bits)]
            rows.add(tuple(sum(bit for a, bit in zip(coefficients, positive_bits)
                               if abs(a) == value) for value in alpha))
        assert rows == large[min(expected_sums)]
        complements = {tuple(m - x for m, x in zip((2, 1, 1, 1, 1), row))
                       for row in rows}
        assert complements == large[max(expected_sums)]
        print(f'Ray {alpha}: five-point fibers {sorted(large)}, '
              f'exact representative u={u}, v={v}, and its complement.')


if __name__ == '__main__':
    main()
