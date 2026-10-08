"""Finite actual-circle checks for linear allocation affine rigidity."""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import atan2, isqrt, pi, sqrt


def mul(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def exact_div(z, w):
    a, b = z
    c, d = w
    norm = c * c + d * d
    x, y = a * c + b * d, b * c - a * d
    if x % norm or y % norm:
        return None
    return x // norm, y // norm


def power(z, n):
    value = (1, 0)
    for _ in range(n):
        value = mul(value, z)
    return value


def factors(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            result.append((p, e))
        p += 1
    if n > 1:
        result.append((n, 1))
    return result


def rank(matrix):
    pivots = {}
    for original in matrix:
        row = list(map(Fraction, original))
        for pivot, previous in sorted(pivots.items()):
            multiple = row[pivot]
            if multiple:
                row = [a - multiple * b for a, b in zip(row, previous)]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            multiple = row[pivot]
            pivots[pivot] = [value / multiple for value in row]
    return len(pivots)


def allocation_data(z, rational_factors):
    allocated, canonical, split = [], (1, 0), []
    for p, e in rational_factors:
        if p == 2:
            canonical = mul(canonical, power((1, 1), e))
        elif p % 4 == 3:
            assert e % 2 == 0
            canonical = mul(canonical, (p ** (e // 2), 0))
        else:
            gaussian = next((a, isqrt(p - a * a)) for a in range(1, isqrt(p) + 1)
                            if isqrt(p - a * a) ** 2 == p - a * a)
            remainder, exponent = z, 0
            while True:
                quotient = exact_div(remainder, gaussian)
                if quotient is None:
                    break
                remainder = quotient
                exponent += 1
            assert 0 <= exponent <= e
            allocated.append(exponent)
            split.append(p)
            canonical = mul(canonical, power(gaussian, exponent))
            canonical = mul(canonical, power((gaussian[0], -gaussian[1]), e - exponent))
    unit = exact_div(z, canonical)
    assert unit in ((1, 0), (-1, 0), (0, 1), (0, -1))
    return tuple(allocated), unit, split


def inspect_arcs(points, allocations, n, constant):
    if not points:
        return 0, 0
    ordered = sorted(points, key=lambda z: atan2(z[1], z[0]) % (2 * pi))
    angles = [atan2(z[1], z[0]) % (2 * pi) for z in ordered]
    doubled = angles + [angle + 2 * pi for angle in angles]
    width = constant / n ** 0.25
    checked, largest = 0, 0
    for start in range(len(ordered)):
        end = start + 1
        # A tiny inward margin avoids treating a floating boundary as
        # evidence for inclusion. These checks do not prove the endpoint.
        while end < start + len(ordered) and doubled[end] - doubled[start] <= width - 1e-12:
            end += 1
        selected = [ordered[j % len(ordered)] for j in range(start, end)]
        rows = [allocations[z] for z in selected]
        assert rank([[1] + list(row) for row in rows]) == len(rows)
        varying = sum(len({row[j] for row in rows}) > 1 for j in range(len(rows[0])))
        assert len(rows) <= varying + 1
        checked += 1
        largest = max(largest, len(rows))
    return checked, largest


def main():
    for e in range(1, 13):
        for a in range(e + 1):
            for b in range(e + 1):
                assert abs(a - b) <= Fraction(a + b) - Fraction(2 * a * b, e)

    maximum_norm = 2500
    circles = defaultdict(list)
    for a in range(-isqrt(maximum_norm), isqrt(maximum_norm) + 1):
        for b in range(-isqrt(maximum_norm), isqrt(maximum_norm) + 1):
            norm = a * a + b * b
            if 1 <= norm <= maximum_norm:
                circles[norm].append((a, b))
    pair_count = arc_count = common_arc_count = 0
    maximum_all = maximum_common = 0
    for n, points in sorted(circles.items()):
        fact = factors(n)
        allocation, unit_groups = {}, defaultdict(list)
        units = {}
        for z in points:
            row, unit, split = allocation_data(z, fact)
            allocation[z] = row
            units[z] = unit
            unit_groups[unit].append(z)
        for z, w in combinations(points, 2):
            pair_norm = 1
            for p, a, b in zip(split, allocation[z], allocation[w]):
                pair_norm *= p ** abs(a - b)
            chord_square = (z[0] - w[0]) ** 2 + (z[1] - w[1]) ** 2
            assert chord_square * pair_norm >= 2 * n
            if units[z] == units[w]:
                assert chord_square * pair_norm >= 4 * n
            pair_count += 1
        count, maximum = inspect_arcs(points, allocation, n, sqrt(2))
        arc_count += count
        maximum_all = max(maximum_all, maximum)
        for group in unit_groups.values():
            count, maximum = inspect_arcs(group, allocation, n, 2)
            common_arc_count += count
            maximum_common = max(maximum_common, maximum)
    print(f"Passed: {pair_count} exact pair tests on {len(circles)} actual circles, N<=2500.")
    print(f"Passed: {arc_count} arbitrary-unit arcs and {common_arc_count} common-unit arcs.")
    print(f"Largest checked arc sizes: {maximum_all} and {maximum_common}; all affine ranks are full.")
    print("All-order product rigidity and the exact boundary constants are proved in the note.")


if __name__ == "__main__":
    main()
