#!/usr/bin/env python3
"""Exact signed common-twist triple identity; no short-arc realization claim."""

from itertools import combinations, product
from fractions import Fraction
from random import Random

PRIMES = ((2, 1, 5), (3, 2, 13), (4, 1, 17))


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def prod(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def primitive(exponents, primes):
    answer = (1, 0)
    for exponent, (x, y, _) in zip(exponents, primes):
        base = (x, y if exponent >= 0 else -y)
        for _ in range(abs(exponent)):
            answer = mul(answer, base)
    return answer


def det3(rows):
    a, b, c = rows
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def determinant(rows):
    rows = [[Fraction(x) for x in row] for row in rows]
    answer = Fraction(1)
    for j in range(len(rows)):
        pivot = next((i for i in range(j, len(rows)) if rows[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            answer = -answer
        value = rows[j][j]
        answer *= value
        for i in range(j + 1, len(rows)):
            scale = rows[i][j] / value
            rows[i] = [x - scale * y for x, y in zip(rows[i], rows[j])]
    assert answer.denominator == 1
    return answer.numerator


def check_general(exponents, primes):
    m = len(exponents)
    d = (m - 1) // 2
    assert m == 2 * d + 1
    parity = [column[0] % 2 for column in zip(*exponents)]
    assert all(all(c % 2 == k for c in column)
               for column, k in zip(zip(*exponents), parity))
    twist = prod(p ** k for (_, _, p), k in zip(primes, parity))
    heights = [prod(p ** ((abs(c) - k) // 2)
                    for c, k, (_, _, p) in zip(row, parity, primes))
               for row in exponents]
    values = [primitive(row, primes) for row in exponents]
    rows = []
    for h, z in zip(heights, values):
        row, power = [h ** d], (1, 0)
        for r in range(1, d + 1):
            power = mul(power, z)
            row.extend([h ** (d - r) * coordinate for coordinate in power])
        rows.append(row)
    pair_values = [primitive([(a - b) // 2 for a, b in zip(exponents[i], exponents[j])], primes)
                   for i, j in combinations(range(m), 2)]
    factor = 1
    for column, k, (_, _, p) in zip(zip(*exponents), parity, primes):
        four_e = (2 * d * sum(map(abs, column))
                  - sum(abs(a - b) for a, b in combinations(column, 2))
                  - 2 * k * d * d)
        same_sign = sum(min(abs(a), abs(b)) for a, b in combinations(column, 2) if a * b > 0)
        assert four_e == 2 * (same_sign - k * d * d)
        assert four_e >= 0 and four_e % 4 == 0
        factor *= p ** (four_e // 4)
    value = determinant(rows)
    assert abs(value) == 2 ** (2 * d * d) * factor * prod(abs(z[1]) for z in pair_values)
    assert factor ** 2 * prod(norm(z) for z in pair_values) == twist ** (d * (d + 1)) * prod(heights) ** (2 * d)
    if d == 1:
        assert value == -4 * factor * prod(z[1] for z in pair_values)


def check(exponents, primes):
    parity = [column[0] % 2 for column in zip(*exponents)]
    assert all(all(c % 2 == k for c in column)
               for column, k in zip(zip(*exponents), parity))
    twist = prod(p ** k for (_, _, p), k in zip(primes, parity))
    heights = [prod(p ** ((abs(c) - k) // 2)
                    for c, k, (_, _, p) in zip(row, parity, primes))
               for row in exponents]
    values = [primitive(row, primes) for row in exponents]
    assert all(norm(z) == twist * h * h for z, h in zip(values, heights))
    pair_values = []
    for r, s in combinations(range(3), 2):
        pair_values.append(primitive([(a - b) // 2
                                     for a, b in zip(exponents[r], exponents[s])], primes))
    factor = 1
    for column, k, (_, _, p) in zip(zip(*exponents), parity, primes):
        numerator = sum(map(abs, column)) - (max(column) - min(column)) - k
        assert numerator >= 0 and numerator % 2 == 0
        factor *= p ** (numerator // 2)
    determinant = det3([(h, z[0], z[1]) for h, z in zip(heights, values)])
    assert determinant == -4 * factor * prod(z[1] for z in pair_values)
    assert factor ** 2 * prod(norm(z) for z in pair_values) == twist ** 2 * prod(heights) ** 2


def main():
    single = 0
    for k in (0, 1):
        candidates = [c for c in range(-7, 8) if c % 2 == k]
        for triple in product(candidates, repeat=3):
            check([(c,) for c in triple], PRIMES[:1])
            single += 1
    rng = Random(20260921)
    for _ in range(100):
        parity = [rng.randrange(2) for _ in PRIMES]
        exponents = [[2 * rng.randrange(-4, 5) + k for k in parity]
                     for _ in range(3)]
        check(exponents, PRIMES)

    characters = [tuple(1 if i in chosen else -1 for i in range(6))
                  for chosen in combinations(range(6), 3) if 0 in chosen]
    source_count = 0
    for _ in range(12):
        allocations = [[rng.randrange(6) for _ in range(6)] for _ in PRIMES]
        data = [tuple(sum(lam[i] * allocation[i] for i in range(6))
                      for allocation in allocations) for lam in characters]
        for triple in combinations(data, 3):
            check(triple, PRIMES)
            source_count += 1
    assert source_count == 1440
    for d in (1, 2):
        for _ in range(40):
            parity = [rng.randrange(2) for _ in PRIMES]
            exponents = [[2 * rng.randrange(-3, 4) + k for k in parity]
                         for _ in range(2 * d + 1)]
            check_general(exponents, PRIMES)
        for positives in range(2 * d + 2):
            count = positives * (positives - 1) // 2
            negatives = 2 * d + 1 - positives
            count += negatives * (negatives - 1) // 2
            assert count == d * d + (positives - d) * (positives - d - 1)
            assert count >= d * d and (count - d * d) % 2 == 0
    unit_fixture = [(1, 3, -1), (-1, 1, 5), (3, -3, 1),
                    (5, -1, -3), (-3, 5, 3)]
    for powers in product(range(4), repeat=len(PRIMES)):
        rotated = []
        for (x, y, p), power in zip(PRIMES, powers):
            for _ in range(power):
                x, y = -y, x
            rotated.append((x, y, p))
        check_general(unit_fixture[:3], rotated)
        check_general(unit_fixture, rotated)
    print('PASS: {} single-prime triples, 100 several-prime triples,'.format(single))
    print('      1440 six-source character triples; exact signed determinant and content.')
    print('PASS: 80 general odd trigonometric determinant/content fixtures at d=1,2.')
    print('PASS: all 64 independent Gaussian-prime unit choices at both d=1 and d=2.')
    print('SCOPE: algebra fixtures only; no short-arc realization or new height exponent.')


if __name__ == '__main__':
    main()
