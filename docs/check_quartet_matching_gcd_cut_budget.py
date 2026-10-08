"""Exact audit of quartet gcds, nested cuts, and primitive fair rows."""

from itertools import combinations
from math import comb, gcd, lcm, prod
from random import Random


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def norm(a):
    return a[0] * a[0] + a[1] * a[1]


def ggcd(a, b):
    while b != (0, 0):
        d = norm(b)
        re = a[0] * b[0] + a[1] * b[1]
        im = a[1] * b[0] - a[0] * b[1]
        quotient = ((re + d // 2) // d, (im + d // 2) // d)
        a, b = b, sub(a, mul(quotient, b))
    return a


def valuation(a):
    """Valuation at 2+i of a nonzero Gaussian integer."""
    assert a != (0, 0)
    count = 0
    while True:
        re, im = 2 * a[0] + a[1], 2 * a[1] - a[0]
        if re % 5 or im % 5:
            return count
        a = re // 5, im // 5
        count += 1


def matching_content(points):
    a, b, c, d = points
    first = mul(sub(b, a), sub(d, c))
    second = mul(sub(c, a), sub(d, b))
    third = mul(sub(d, a), sub(c, b))
    assert sub(second, first) == third
    return ggcd(first, second)


def audit_global_gcd():
    rng = Random(9142026)
    pool = [(x, y) for x in range(-12, 13) for y in range(-12, 13)]
    for _ in range(1000):
        points = rng.sample(pool, rng.randrange(4, 9))
        delta = {(i, j): valuation(sub(points[j], points[i]))
                 for i, j in combinations(range(len(points)), 2)}
        r = min(delta.values())
        # Same class iff the pair valuation exceeds r.
        classes = []
        for i in range(len(points)):
            for group in classes:
                if delta[tuple(sorted((i, group[0])))] > r:
                    group.append(i)
                    break
            else:
                classes.append([i])
        largest = max(classes, key=len)
        expected = 2 * r
        if len(largest) == len(points) - 1:
            expected = r + min(delta[tuple(sorted(pair))]
                               for pair in combinations(largest, 2))
        content = (0, 0)
        for quartet in combinations(points, 4):
            content = ggcd(content, matching_content(quartet))
        assert valuation(content) == expected


def binom(n, k):
    return comb(n, k) if n >= k >= 0 else 0


def J(m, s):
    return (2 * (binom(s, 4) + binom(m - s, 4))
            + (m - s) * binom(s, 3) + s * binom(m - s, 3))


def audit_cuts():
    for m in range(4, 17):
        for s in range(m + 1):
            direct = sum(abs(len(set(q) & set(range(s))) - 2)
                         for q in combinations(range(m), 4))
            assert direct == J(m, s)
            assert J(m, s) >= J(m, m // 2)
            assert (m - 3) * (binom(s, 3) + binom(m - s, 3)) == (
                J(m, s) + 2 * (binom(s, 4) + binom(m - s, 4)))
            if s < m:
                t = m - s - 1
                assert 6 * (J(m, s + 1) - J(m, s)) == (
                    (s - t) * ((m - 2) * (m - 3) + 2 * s * t))
        if m % 2 == 0:
            n = m // 2
            assert J(m, n) * (2 * n - 1) * (2 * n - 3) == (
                comb(m, 4) * 3 * (n - 1) * (n - 2))
    rng = Random(914)
    for _ in range(1000):
        e = rng.randrange(1, 15)
        levels = sorted(rng.randrange(e + 1) for _ in range(4))
        direct = 2 * e + sum(levels[:2]) - sum(levels[2:])
        layers = sum(abs(sum(t >= h for t in levels) - 2)
                     for h in range(1, e + 1))
        assert direct == layers


def audit_actual_rows():
    m = 6
    cuts = [set(cut) for cut in combinations(range(m), m // 2)]
    K = len(cuts)
    M = lcm(*(j * j - k * k for j in range(1, K + 1) for k in range(1, j)))
    for q in (1, 2):
        blocks = [(2 * j * M * q, 1) for j in range(1, K + 1)]
        norms = [norm(block) for block in blocks]
        assert all(gcd(a, b) == 1 for a, b in combinations(norms, 2))
        points = []
        for i in range(m):
            z = (1, 0)
            for block, cut in zip(blocks, cuts):
                z = mul(z, block if i in cut else (block[0], -block[1]))
            points.append(z)
        N = prod(norms)
        assert len(set(points)) == m and all(norm(z) == N for z in points)
        row_gcd = (0, 0)
        for z in points:
            row_gcd = ggcd(row_gcd, z)
        assert norm(row_gcd) == 1
        all_gamma = (0, 0)
        for quartet in combinations(range(m), 4):
            gamma = matching_content([points[i] for i in quartet])
            forced = prod(n ** abs(len(set(quartet) & cut) - 2)
                          for n, cut in zip(norms, cuts))
            assert norm(gamma) % forced == 0
            all_gamma = ggcd(all_gamma, gamma)
        assert gcd(norm(all_gamma), N) == 1


if __name__ == "__main__":
    audit_global_gcd()
    audit_cuts()
    audit_actual_rows()
    print("PASS: 1,000 global gcd cases, cut identities through m=16,")
    print("1,000 nested profiles, and two exact primitive fair six-row realizations.")
