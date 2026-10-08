"""Finite audit of the weighted four-cube five-point gap."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import atan2, pi, sqrt
from random import Random


def det(matrix):
    if not matrix:
        return 1
    return sum((-1) ** j * value * det([
        row[:j] + row[j + 1:] for row in matrix[1:]
    ]) for j, value in enumerate(matrix[0]))


def adj(matrix):
    r = len(matrix)
    return [[(-1) ** (i + j) * det([
        row[:i] + row[i + 1:] for k, row in enumerate(matrix) if k != j
    ]) for j in range(r)] for i in range(r)]


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def product_int(values):
    result = 1
    for value in values:
        result *= value
    return result


def gaussian_audit(rows, unit_exponents, common):
    gaussian = [(2, 1), (3, 2), (4, 1), (5, 2)]
    primes = list(map(norm, gaussian))
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    points = []
    for row, unit in zip(rows, unit_exponents):
        z = mul(common, units[unit])
        for bit, (a, b) in zip(row, gaussian):
            z = mul(z, (a, b if bit else -b))
        points.append(z)
    full_common = common
    varying = []
    for j, (a, b) in enumerate(gaussian):
        values = {row[j] for row in rows}
        if len(values) == 1:
            full_common = mul(full_common, (a, b if rows[0][j] else -b))
        else:
            varying.append(primes[j])
    n0 = product_int(varying)
    n = norm(full_common) * n0
    assert all(norm(z) == n for z in points)
    good_pairs = []
    for i, j in combinations(range(5), 2):
        p = product_int(prime for prime, x, y in zip(primes, rows[i], rows[j]) if x != y)
        if p * p <= n0:
            good_pairs.append((i, j))
            chord2 = norm((points[i][0] - points[j][0], points[i][1] - points[j][1]))
            constant = 4 if unit_exponents[i] == unit_exponents[j] else 2
            assert chord2 * p >= constant * n
            # Squared version of chord2 >= 2 R |d|, avoiding radicals.
            assert chord2 ** 2 >= 4 * n * norm(full_common)
            if len(set(unit_exponents)) == 1:
                assert chord2 ** 2 >= 16 * n * norm(full_common)
    assert good_pairs
    return len(good_pairs)


def phase_bound_audit():
    targets = [pi * m / (4 * h) for h in [1, 2, 3] for m in range(8 * h)]
    checked = 0
    for a in range(1, 60):
        for b in range(1, 60):
            p = a * a + b * b
            if p % 4 != 1 or any(p % d == 0 for d in range(2, int(sqrt(p)) + 1)):
                continue
            phi = atan2(b, a)
            distance = min(abs((phi - beta + pi) % (2 * pi) - pi) for beta in targets)
            assert distance >= 1 / (4 * p) - 1e-14
            checked += 1
    return checked


def main():
    vertices = list(product([0, 1], repeat=4))
    weights = [(1, 1, 1, 1), (1, 2, 3, 5), (17, 1, 1, 1),
               (0, 1, 2, 3), (0, 0, 0, 0), (997, 1009, 1013, 1021)]
    counts = Counter()
    survivors = []
    all_rows = list(combinations(vertices, 5))
    max_cost = Fraction(0)
    for rows in all_rows:
        distances = [hamming(a, b) for a, b in combinations(rows, 2)]
        if 1 not in distances and 4 not in distances:
            centers = [row for row in rows if all(other == row or hamming(row, other) == 3
                                                  for other in rows)]
            assert len(centers) == 1
            leaves = [row for row in rows if row != centers[0]]
            cut_counts = [sum(a[j] != b[j] for a, b in combinations(leaves, 2))
                          for j in range(4)]
            assert cut_counts == [3, 3, 3, 3]
            survivors.append(rows)
        for weight in weights:
            assert any(2 * sum(w for w, x, y in zip(weight, a, b) if x != y) <= sum(weight)
                       for a, b in combinations(rows, 2))
        matrix = [[a - b for a, b in zip(row, rows[0])] for row in rows[1:]]
        delta = det(matrix)
        counts[abs(delta)] += 1
        if delta:
            assert abs(delta) <= 3
            jmat = adj(matrix)
            assert all(abs(value) <= 2 for row in jmat for value in row)
            for row in jmat:
                # Exact best error bound from the source containing interval.
                cost = Fraction(sum(map(abs, row)) + abs(sum(row)), 4 * abs(delta))
                max_cost = max(max_cost, cost)
                assert cost <= Fraction(3, 2)
    assert counts == {0: 1360, 1: 2672, 2: 320, 3: 16}
    assert len(survivors) == 16
    assert max_cost == Fraction(3, 2)
    rng = Random(20260915)
    actual = 0
    for rows in survivors + rng.sample(all_rows, 64):
        actual += gaussian_audit(rows, [rng.randrange(4) for _ in rows], (7, 11))
        actual += gaussian_audit(rows, [0] * 5, (2, 3))
    phases = phase_bound_audit()
    print(f"PASS: {len(all_rows)} five-vertex subsets; {len(survivors)} translated "
          f"four-triple configurations; determinant counts {dict(sorted(counts.items()))}.")
    print(f"Maximum exact inverse phase cost {max_cost}; {actual} literal Gaussian "
          f"pair certificates; {phases} quadratic-grid prime checks.")


if __name__ == '__main__':
    main()
