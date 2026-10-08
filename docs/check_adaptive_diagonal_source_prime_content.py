"""Exact source-prime clipping checks for all source-anchor target radii."""

from itertools import combinations
from math import gcd
from random import Random

from check_adaptive_diagonal_anchor_radius_product import (
    check, conj, det, ggcd, mul, norm, split_prime, vg, vp,
)


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def clipping_fixture(rows, a, d):
    assert a != d
    source_n, big_b, edges, target_n, quotient = check(rows, a, d)
    m, u, v = len(rows), a + d, a - d
    overlap = gcd(source_n, abs(u * v))
    assert big_b ** (m - 2) * overlap ** ((m - 1) * (m - 2)) % quotient == 0
    beta = {(i, k): abs(det(rows[i], rows[k])) // norm(ggcd(rows[i], rows[k]))
            for i, k in combinations(range(m), 2)}
    local_checks = spikes = nonzero_levels = 0

    for p in (5, 13, 17, 29, 37):
        pi = split_prime(p)
        source = [vg(z, pi) - vg(z, conj(pi)) for z in rows]
        width = max(source) - min(source)
        assert width == vp(source_n, p)
        h = vp(abs(u * v), p)
        baseline_loss = actual_loss = 0
        for j in range(m):
            actual, baseline, excess = [], [], []
            for i in range(m):
                w = mul(rows[i], conj(rows[j]))
                target = a * w[0], d * w[1]
                t = vg(target, pi) - vg(target, conj(pi))
                x = source[i] - source[j]
                actual.append(t)
                if h == 0:
                    if x:
                        assert t == 0
                    baseline.append(0)
                    excess.append(0)
                else:
                    sign = 1 if v % p == 0 else -1
                    clipped = sign * max(-h, min(x, h))
                    baseline.append(clipped)
                    if abs(x) != h:
                        assert t == clipped
                        excess.append(0)
                    else:
                        assert t * clipped > 0 and abs(t) >= h
                        excess.append(abs(t) - h)
                        spikes += abs(t) > h
                nonzero_levels += bool(width and t)

            row_loss = sum(map(abs, actual)) - max(actual) + min(actual)
            base_loss = sum(map(abs, baseline)) - max(baseline) + min(baseline)
            assert row_loss >= base_loss
            assert vp(target_n[j], p) == max(actual) - min(actual)
            local_budget = 0
            for i, k in combinations([i for i in range(m) if i != j], 2):
                budget = vp(beta[min(i, k), max(i, k)], p)
                local_budget += budget
                if actual[i] * actual[k] > 0:
                    if h == 0:
                        assert source[i] == source[k] == source[j]
                        assert min(abs(actual[i]), abs(actual[k])) <= budget
                    elif excess[i] and excess[k]:
                        assert source[i] == source[k]
                        assert min(excess[i], excess[k]) <= budget
            assert row_loss <= base_loss + local_budget
            baseline_loss += base_loss
            actual_loss += row_loss

        predicted = (2 * sum(min(abs(source[i] - source[j]), h)
                             for i, j in combinations(range(m), 2))
                     - sum(min(s - min(source), h) + min(max(source) - s, h)
                           for s in source))
        assert predicted == baseline_loss >= 0
        if width:
            low, high = source.index(min(source)), source.index(max(source))
            assert baseline_loss == sum(
                min(abs(source[i] - source[j]), h)
                for i in range(m) if i not in (low, high) for j in range(m))
        assert baseline_loss <= (m - 1) * (m - 2) * min(h, width)
        assert actual_loss == vp(quotient, p)
        assert actual_loss <= baseline_loss + (m - 2) * vp(big_b, p)
        local_checks += 1
    return local_checks, spikes, nonzero_levels, bool(overlap > 1), bool(quotient > 1)


def main():
    fixture = check([(1, 0), (2, 1), (3, 1)], 3, 2)
    assert fixture == (25, 1, [5, 85, 445], [425, 445, 7565], 25)
    assert fixture[-1] > fixture[1]  # The overlap correction is necessary.

    rows_fixtures = [
        [(1, 0), (2, 1), (3, 1)],
        [(1, 0), (1, 1), (2, 1), (1, 2)],
        [(1, 0), (4, 1), (9, 1), (14, 1), (19, 1)],
    ]
    for pi in ((2, 1), (3, 2)):
        for levels in ((0, 1, -1), (0, 1, 2, 3), (0, -2, -1, 1, 2, 4)):
            rows_fixtures.append([power(pi if s >= 0 else conj(pi), abs(s))
                                  for s in levels])
    rng = Random(340)
    for _ in range(60):
        rows = [(1, 0)]
        wanted = rng.randrange(3, 7)
        while len(rows) < wanted:
            z = rng.randrange(-15, 16), rng.randrange(1, 16)
            if gcd(*z) == 1 and all(det(z, w) for w in rows):
                rows.append(z)
        rows_fixtures.append(rows)

    maps = ((2, 1), (1, 2), (3, 1), (3, 2), (6, 1), (13, 12),
            (26, 1), (7, 6), (14, 1), (85, 84), (170, 1),
            (5, 3), (9, 7), (17, 15))
    totals = [0] * 5
    count = 0
    for rows in rows_fixtures:
        for a, d in maps:
            result = clipping_fixture(rows, a, d)
            totals = [x + y for x, y in zip(totals, result)]
            if (a, d) in ((2, 1), (3, 1), (5, 3), (9, 7), (17, 15)):
                n, b, _, _, quotient = check(rows, a, d)
                assert gcd(n, (a + d) * (a - d)) == 1
                assert b ** (len(rows) - 2) % quotient == 0
            count += 1
    assert totals[1] > 0 and totals[3] > 0
    print("PASS: %d global source-prime divisibilities; %d local clipping audits; "
          "%d outward boundary spikes; %d nonzero target valuations at source primes; "
          "%d source-parameter overlaps; %d nontrivial quotients."
          % (count, *totals))


if __name__ == "__main__":
    main()
