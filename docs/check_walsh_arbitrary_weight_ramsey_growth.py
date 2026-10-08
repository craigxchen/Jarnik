"""Bounded exact identities behind the arbitrary-weight Ramsey growth proof."""

from fractions import Fraction as Q
from itertools import combinations
from random import Random


def dot(a, b):
    return bin(a & b).count("1") % 2


def averages():
    rng = Random(333)
    count = 0
    for m in (8, 16):
        labels = list(range(1, m)) + [1]
        rng.shuffle(labels)
        weights = [Q(rng.randrange(1, 31), 7) for _ in range(m)]
        total = sum(weights)
        planes = {tuple(sorted((0, d, e, d ^ e)))
                  for d, e in combinations(range(1, m), 2)}
        for plane in planes:
            _, d, e, _ = plane
            assert {0, d, e, d ^ e} == set(plane)
            seen = set()
            cosets = []
            for x in range(m):
                if x not in seen:
                    rows = {x ^ v for v in plane}
                    cosets.append(rows)
                    seen |= rows
            for ell in (1, 2, 3):
                matches = [dot(a, d) + 2 * dot(a, e) == ell for a in labels]
                cell = sum(w for w, good in zip(weights, matches) if good)
                values = []
                for rows in cosets:
                    fp = sum(weights[x] for x in rows)
                    fm = sum(weights[x] for x in rows if matches[x])
                    values.append(2 * (2 * fm - fp) + Q(4, m) * (total - 4 * cell))
                # The separate analytic term is -4 log 2 in every coset.
                assert sum(values) / len(values) == -Q(4, m) * total
                count += 1
    return count


def slack_and_parseval():
    rng = Random(334)
    triangles = spectra = 0
    for m in (4, 8, 16, 32):
        for _ in range(20):
            u = [Q(0)] + [Q(rng.randrange(1, 20), 5) for _ in range(m - 1)]
            total, error = sum(u), Q(8)
            for d in range(1, m):
                margins = [u[e] + u[d ^ e] + error - u[d]
                           for e in range(1, m) if e != d]
                assert min(margins) >= 0
                assert sum(margins) == 2 * total + (m - 2) * error - m * u[d]
                assert u[d] <= 2 * total / m + Q(m - 2, m) * error
                triangles += 1

            q = [Q(0)] + [Q(rng.randrange(1, 40), 7) for _ in range(m - 1)]
            qsum = sum(q)
            transform = [sum(q[a] * (-1) ** dot(a, d) for a in range(m))
                         for d in range(m)]
            kappa = max(transform[1:]) + Q(1)
            slack = [kappa - value for value in transform[1:]]
            mean = sum(slack) / (m - 1)
            centered = [Q(0)] + [x - qsum / (m - 1) for x in q[1:]]
            assert sum(slack) == qsum + (m - 1) * kappa
            assert all(sum(centered[a] * (-1) ** dot(a, d) for a in range(m))
                       == mean - slack[d - 1] for d in range(1, m))
            variance_sum = sum((x - mean) ** 2 for x in slack)
            assert sum(x * x for x in centered) == variance_sum / m
            bound = max(slack)
            assert variance_sum <= Q(m - 1, 4) * bound * bound
            assert sum(map(abs, centered)) ** 2 <= Q((m - 1) ** 2, 4 * m) * bound ** 2
            spectra += 1
    return triangles, spectra


def small_ramsey():
    planes = {tuple(sorted((d, e, d ^ e)))
              for d, e in combinations(range(1, 8), 2)}
    assert len(planes) == 7
    for colors in range(128):
        assert any(len({(colors >> (d - 1)) & 1 for d in plane}) == 1 for plane in planes)
    # The general dimensional claim is an external theorem, not this finite test.
    return 128


def tighter_parameters():
    # Integer bounds suffice: 15 + 2 log_2(m) <= 2m follows from
    # 2^15 m^2 <= 2^(2m), so h >= s/(144m^2) > 32K.
    tested = 0
    for exponent in range(24, 257):
        for k in (2 ** exponent, 2 ** exponent + 1,
                  2 ** (exponent + 1) - 1):
            m = (k - 1).bit_length()
            s = 2 ** 14 * k * m * m
            assert 2 ** 15 * m * m <= 2 ** (2 * m)
            assert s <= k * k
            assert Q(s, 144 * m * m) > 32 * k
            assert 50 * 2 ** 14 < 2 ** 20
            tested += 1

    # Test T=2^j without evaluating any astronomical tower in M.
    # The integer square root makes K's floor exact also for odd j.
    from math import isqrt
    for j in range(100, 513):
        tstar = 2 ** j
        k = isqrt(tstar) // (2 ** 12 * j)
        assert k >= 2 ** 24
        m = (k - 1).bit_length()
        s = 2 ** 14 * k * m * m
        assert m <= j
        assert s <= k * k <= tstar
        assert 2 ** 20 * k * k * m * m <= tstar // 16
        assert tstar - 1 > tstar // 16
        tested += 1
    return tested


if __name__ == "__main__":
    row_averages = averages()
    triangles, spectra = slack_and_parseval()
    colorings = small_ramsey()
    parameters = tighter_parameters()
    assert 2 ** 24 >= 512 * (2 * 24 + 1) ** 2
    print("PASS: %d exact four-row averages; %d slack-maximum identities; "
          "%d Fourier/variance fixtures; %d finite two-colorings; "
          "%d refined parameter fixtures."
          % (row_averages, triangles, spectra, colorings, parameters))
