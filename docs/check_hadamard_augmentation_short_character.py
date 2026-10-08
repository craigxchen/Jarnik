"""Exact short characters for complete Hadamard cores with extra columns.

This checks combinatorics and certificates, not endpoint realizations or
the linear-form theorem used in the companion effective phase note.
"""

from itertools import combinations, product
from math import comb, isqrt
from random import Random

from check_rank_one_kernel_relative_roth import paley_core_11


def walsh(order):
    assert order >= 4 and order & (order - 1) == 0
    return [[-1 if bin(i & j).count("1") % 2 else 1
             for j in range(order)] for i in range(order)]


def paley():
    core = paley_core_11()
    return [[1] * 12] + [[1] + row for row in core]


def canonical(column):
    return tuple(value * column[0] for value in column)


def audit_hadamard(h):
    m = len(h)
    assert all(len(row) == m and row[0] == 1 for row in h)
    assert all(value in (-1, 1) for row in h for value in row)
    assert all(sum(h[i][j] * h[i][l] for i in range(m)) == m * (j == l)
               for j in range(m) for l in range(m))


def fourier(h, extras):
    m, k = len(h), len(extras)
    numerators = [[sum(h[i][j] * extras[l][i] for i in range(m))
                   for l in range(k)] for j in range(1, m)]
    means = [sum(column) for column in extras]
    for l in range(k):
        assert means[l] ** 2 + sum(row[l] ** 2 for row in numerators) == m * m
    assert sum(x * x for row in numerators for x in row) <= k * m * m
    good = [j for j, row in enumerate(numerators)
            if (m - 1) * sum(x * x for x in row) <= 2 * k * m * m]
    assert len(good) >= m // 2
    bound = isqrt(2 * k * m * m // (m - 1))
    assert all(abs(x) <= bound for j in good for x in numerators[j])
    return numerators, good


def verify_certificate(h, extras, v):
    m = len(h)
    assert len(v) == m - 1 and any(v)
    assert set(v) <= {-1, 0, 1}
    lam = [sum(h[i][j + 1] * v[j] for j in range(m - 1)) for i in range(m)]
    ell = sum(map(abs, v))
    assert sum(lam) == 0
    assert sum(x * x for x in lam) == m * ell
    assert sum(map(abs, lam)) ** 2 <= m * m * ell
    columns = list(zip(*h))[1:] + list(extras)
    signature = [sum(lam[i] * column[i] for i in range(m)) for column in columns]
    assert signature == [m * x for x in v] + [0] * len(extras)
    # The signed Gaussian exponents are half this signature. Common
    # content cancels because sum(lam)=0, with units still external.
    assert [sum(lam[i] * ((column[i] + 1) // 2) for i in range(m))
            for column in columns] == [m * x // 2 for x in v] + [0] * len(extras)
    return ell


def one_extra_certificate(numerators):
    seen = {}
    for j, (value,) in enumerate(numerators):
        v = [0] * len(numerators)
        if value == 0:
            v[j] = 1
            return v
        if abs(value) in seen:
            prior, prior_value = seen[abs(value)]
            v[j], v[prior] = 1, -value // prior_value
            return v
        seen[abs(value)] = (j, value)
    raise AssertionError("Parseval forces a zero or equal absolute numerators")


def subset_collision(numerators, candidates, q):
    seen = {}
    k = len(numerators[0])
    for subset in combinations(candidates, q):
        value = tuple(sum(numerators[j][l] for j in subset) for l in range(k))
        if value in seen:
            v = [0] * len(numerators)
            for j in subset:
                v[j] += 1
            for j in seen[value]:
                v[j] -= 1
            assert any(v) and sum(map(abs, v)) <= 2 * q
            return v
        seen[value] = subset
    return None


def exhaustive_one_extra():
    checked = 0
    h = walsh(8)
    audit_hadamard(h)
    old = {canonical(column) for column in zip(*h)}
    for tail in product((-1, 1), repeat=7):
        extra = (1,) + tail
        if canonical(extra) in old:
            continue
        numerators, _ = fourier(h, [extra])
        v = one_extra_certificate(numerators)
        assert verify_certificate(h, [extra], v) <= 2
        checked += 1
    assert checked == 120
    return checked


def random_fixtures():
    rng = Random(261)
    checked = 0
    for h in (walsh(16), walsh(32), walsh(64), paley()):
        audit_hadamard(h)
        m = len(h)
        for k in (1, 2, 3, 4):
            old = {canonical(column) for column in zip(*h)}
            extras = []
            while len(extras) < k:
                column = tuple([1] + [rng.choice((-1, 1)) for _ in range(m - 1)])
                if canonical(column) not in old:
                    old.add(canonical(column))
                    extras.append(column)
            numerators, good = fourier(h, extras)
            if k == 1:
                v = one_extra_certificate(numerators)
            else:
                v = None
                for q in range(1, 5):
                    v = subset_collision(numerators, good, q)
                    if v is not None:
                        break
            assert v is not None, (m, k)
            verify_certificate(h, extras, v)
            # These profiles pass the abstract strict Gram test with old
            # weights 1 and extra weights 1/(2k).
            assert all(2 * k * sum(h[i][j] * h[l][j] for j in range(1, m))
                       + sum(f[i] * f[l] for f in extras) < 0
                       for i, l in combinations(range(m), 2))
            checked += 1
    return checked


def pigeonhole_bounds():
    checked = 0
    for k in range(1, 65):
        a = (2 * k - 1).bit_length()
        assert 2 ** (a - 1) < 2 * k <= 2 ** a
        q = 4 * k * a
        assert 2 ** (24 * a) > 160 ** 2 * k ** 4 * a ** 3
        for multiplier in (16, 20, 32):
            m = multiplier * q
            b = isqrt(2 * k * m * m // (m - 1))
            assert comb(m // 2, q) > (2 * q * b + 1) ** k
            checked += 1
        # An explicit coarse threshold for the fixed-k smaller q comes
        # from squaring (M/(2q))^q > (5q sqrt(kM))^k.
        small_q = k // 2 + 1
        exponent = 2 * small_q - k
        threshold_power = ((2 * small_q) ** (2 * small_q)
                           * (5 * small_q) ** (2 * k) * k ** k)
        m = 2 * (threshold_power + 1 if exponent == 1
                 else isqrt(threshold_power) + 1)
        assert exponent in (1, 2) and m >= 2 * small_q
        assert m ** exponent > threshold_power
        b = isqrt(2 * k * m * m // (m - 1))
        assert comb(m // 2, small_q) > (2 * small_q * b + 1) ** k
        checked += 1
    return checked


def truncated_core_obstruction():
    checked = 0
    for n in (8, 16, 32, 64):
        h = walsh(n)
        rows = [row[1:] for row in h[2:]]
        m, t = n - 2, n - 1
        assert t == m + 1 and m % 4 == 2 and m > 2
        columns = list(zip(*rows))
        assert all(len(set(column)) == 2 for column in columns)
        assert len(set(map(canonical, columns))) == t
        assert all(sum(x * y for x, y in zip(rows[i], rows[j])) == -1
                   for i, j in combinations(range(m), 2))
        checked += 1
    return checked


if __name__ == "__main__":
    exhaustive = exhaustive_one_extra()
    fixtures = random_fixtures()
    bounds = pigeonhole_bounds()
    truncated = truncated_core_obstruction()
    print(f"PASS: {exhaustive} exhaustive one-extra cases and {fixtures} Walsh/Paley fixtures.")
    print("PASS: exact full-matrix integer certificates, Parseval, parity, and signed exponent cancellation.")
    print(f"PASS: {bounds} finite pigeonhole inequalities and {truncated} truncated-core obstructions.")
    print("No actual endpoint angles or general low-redundancy extraction are asserted.")
