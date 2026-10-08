#!/usr/bin/env python3
"""Exact finite checks for exact_nested_profile_residual_extraction.md.

No floating point arithmetic or external packages are used. These tests
verify local prime-layer identities and selected arithmetic consequences;
they do not prove the endpoint tuple-selection theorem.
"""

from __future__ import annotations

from itertools import combinations, product
from math import lcm, prod

Gauss = tuple[int, int]
ZERO: Gauss = (0, 0)
ONE: Gauss = (1, 0)


def mul(a: Gauss, b: Gauss) -> Gauss:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def conj(a: Gauss) -> Gauss:
    return (a[0], -a[1])


def norm(a: Gauss) -> int:
    return a[0] * a[0] + a[1] * a[1]


def power(a: Gauss, e: int) -> Gauss:
    out = ONE
    for _ in range(e):
        out = mul(out, a)
    return out


def nearest(n: int, d: int) -> int:
    q, r = divmod(n, d)
    return q + (2 * r > d)


def rem(a: Gauss, b: Gauss) -> Gauss:
    q0, q1 = mul(a, conj(b))
    den = norm(b)
    q = (nearest(q0, den), nearest(q1, den))
    prod = mul(q, b)
    out = (a[0] - prod[0], a[1] - prod[1])
    assert norm(out) < den
    return out


def gaussian_gcd(a: Gauss, b: Gauss) -> Gauss:
    while b != ZERO:
        a, b = b, rem(a, b)
    return a


def gcd_all(values: list[Gauss]) -> Gauss:
    out = ZERO
    for value in values:
        out = gaussian_gcd(out, value)
    return out


def norm_gcd(values: list[Gauss]) -> int:
    return norm(gcd_all(values))


def build(
    profiles: list[tuple[Gauss, int, tuple[int, ...]]],
) -> tuple[list[Gauss], dict[tuple[int, ...], Gauss], int]:
    """Rows and exact threshold blocks for one anchor plus m other rows."""
    m = len(profiles[0][2]) - 1
    rows = [ONE for _ in range(m)]
    blocks: dict[tuple[int, ...], Gauss] = {}
    empty_mass = 1
    for pi, e, allocations in profiles:
        assert len(allocations) == m + 1
        assert 0 <= min(allocations) <= max(allocations) <= e
        a0 = allocations[0]
        for i, ai in enumerate(allocations[1:]):
            factor = mul(power(pi, max(ai - a0, 0)),
                         power(conj(pi), max(a0 - ai, 0)))
            rows[i] = mul(rows[i], factor)
        signs_by_t: dict[tuple[int, ...], set[int]] = {}
        for t in range(1, e + 1):
            if t <= a0:
                support = tuple(i for i, ai in enumerate(allocations[1:])
                                if ai < t)
                orient = -1
            else:
                support = tuple(i for i, ai in enumerate(allocations[1:])
                                if ai >= t)
                orient = 1
            if not support:
                empty_mass *= norm(pi)
                continue
            signs_by_t.setdefault(support, set()).add(orient)
            blocks[support] = mul(
                blocks.get(support, ONE), pi if orient == 1 else conj(pi)
            )
        assert all(len(orientations) == 1
                   for orientations in signs_by_t.values())
    return rows, blocks, empty_mass


def circle_rows(
    profiles: list[tuple[Gauss, int, tuple[int, ...]]],
) -> list[Gauss]:
    count = len(profiles[0][2])
    rows = [ONE for _ in range(count)]
    for pi, e, allocations in profiles:
        assert len(allocations) == count
        for i, ai in enumerate(allocations):
            rows[i] = mul(
                rows[i], mul(power(pi, ai), power(conj(pi), e - ai))
            )
    expected_norm = prod(norm(pi) ** e for pi, e, _ in profiles)
    assert all(norm(row) == expected_norm for row in rows)
    return rows


def direct_pair(
    profiles: list[tuple[Gauss, int, tuple[int, ...]]], i: int, j: int
) -> Gauss:
    out = ONE
    for pi, _, allocations in profiles:
        ai, aj = allocations[i], allocations[j]
        out = mul(
            out, mul(power(pi, max(ai - aj, 0)),
                     power(conj(pi), max(aj - ai, 0)))
        )
    return out


def check_profile(
    profiles: list[tuple[Gauss, int, tuple[int, ...]]],
) -> tuple[int, int, int]:
    rows, blocks, empty_mass = build(profiles)
    m = len(rows)
    assert norm_gcd(circle_rows(profiles)) == empty_mass
    for i in range(m):
        prod = ONE
        for support, block in blocks.items():
            if i in support:
                prod = mul(prod, block)
        assert prod == rows[i]
        assert norm_gcd([rows[i], conj(rows[i])]) == 1
        assert rows[i] == direct_pair(profiles, i + 1, 0)
    pair_checks = 0
    for i, j in combinations(range(m), 2):
        g = norm_gcd([rows[i], rows[j]])
        a_ij = direct_pair(profiles, i + 1, j + 1)
        left = mul(rows[i], conj(rows[j]))
        assert left == (g * a_ij[0], g * a_ij[1])
        pair_checks += 1
    for r in range(1, m + 1):
        for inds in combinations(range(m), r):
            expected = 1
            for support, block in blocks.items():
                if all(i in support for i in inds):
                    expected *= norm(block)
            assert norm_gcd([rows[i] for i in inds]) == expected
    return 2**m - 1, 1, pair_checks


def check_anchor_invariance(
    profiles: list[tuple[Gauss, int, tuple[int, ...]]],
) -> tuple[int, int]:
    count = len(profiles[0][2])
    expected_values = [
        abs(direct_pair(profiles, i, j)[1])
        for i, j in combinations(range(count), 2)
    ]
    assert all(value > 0 for value in expected_values)
    expected_product = prod(expected_values)
    original_content = norm_gcd(circle_rows(profiles))
    anchor_checks = 0
    residual_checks = 0
    for anchor in range(count):
        order = [anchor] + [i for i in range(count) if i != anchor]
        reanchored = [
            (pi, e, tuple(allocations[i] for i in order))
            for pi, e, allocations in profiles
        ]
        rows, _, empty_mass = build(reanchored)
        assert empty_mass == original_content
        actual_values = [abs(row[1]) for row in rows]
        for i, j in combinations(range(count - 1), 2):
            g = norm_gcd([rows[i], rows[j]])
            delta = mul(conj(rows[i]), rows[j])[1]
            assert delta % g == 0
            actual_values.append(abs(delta // g))
        assert all(value > 0 for value in actual_values)
        assert sorted(actual_values) == sorted(expected_values)
        assert prod(actual_values) == expected_product
        anchor_checks += 1
        residual_checks += len(actual_values)
    return anchor_checks, residual_checks


def check_first_moment() -> int:
    count = 0
    for size in range(2, 6):
        for e in range(1, 5):
            for aa in product(range(e + 1), repeat=size):
                pair_sum = sum(abs(aa[i] - aa[j])
                               for i in range(size) for j in range(i + 1, size))
                layer_sum = sum(
                    sum(a >= t for a in aa) * sum(a < t for a in aa)
                    for t in range(1, e + 1)
                )
                assert pair_sum == layer_sum
                assert 4 * layer_sum <= e * size * size
                count += 1
    return count


def valuation(n: int, p: int) -> int:
    count = 0
    while n % p == 0:
        n //= p
        count += 1
    return count


def hankel_sum(q: list[int], c: list[int], subset_size: int) -> int:
    total = 0
    for inds in combinations(range(len(q)), subset_size):
        selected = set(inds)
        term = 1
        for i, j in combinations(inds, 2):
            term *= (q[i] - q[j]) ** 2
        for i in range(len(q)):
            if i not in selected:
                term *= c[i]
        total += term
    return total


def check_hankel_overlap() -> tuple[int, int]:
    # One rational prime has two disjoint nested chains:
    # {0} subset {0,1} below the anchor and
    # {3} subset {2,3} above it.
    pi = (2, 1)  # Norm 5.
    rows, blocks, empty_mass = build([(pi, 4, (2, 0, 1, 3, 4))])
    assert empty_mass == 1
    assert {s: norm(v) for s, v in blocks.items()} == {
        (0,): 5, (0, 1): 5, (2, 3): 5, (3,): 5
    }
    ts: list[int] = []
    for i, j in combinations(range(4), 2):
        delta = rows[i][0] * rows[j][1] - rows[j][0] * rows[i][1]
        g = norm_gcd([rows[i], rows[j]])
        assert delta and delta % g == 0
        ts.append(abs(delta // g))
    scale = lcm(*(abs(row[1]) for row in rows), *ts)
    q = [scale * x // y for x, y in rows]
    assert all(scale * x % y == 0 for x, y in rows)
    assert len(set(q)) == 4
    c = [x * x + scale * scale for x in q]
    assert all(c[i] == (scale // rows[i][1]) ** 2 * norm(rows[i])
               for i in range(4))

    strict = 0
    checks = 0
    for r in range(1, 5):
        determinant = hankel_sum(q, c, r)
        v = valuation(determinant, 5)
        per_subset = []
        for inds in combinations(range(4), r):
            selected = set(inds)
            val = 0
            for support in blocks:
                t = len(selected.intersection(support))
                val += len(support) - t + t * (t - 1)
            per_subset.append(val)
        full_min = min(per_subset)
        separate_min = sum(
            min(
                len(support) - len(set(inds).intersection(support))
                + len(set(inds).intersection(support))
                * (len(set(inds).intersection(support)) - 1)
                for inds in combinations(range(4), r)
            )
            for support in blocks
        )
        assert v >= full_min >= separate_min
        strict += full_min > separate_min
        checks += 1
    assert strict >= 1
    return checks, strict


def main() -> None:
    local_profiles = 0
    gcd_checks = 0
    original_gcd_checks = 0
    pair_checks = 0
    for m in range(1, 5):
        for e in range(1, 5):
            for allocations in product(range(e + 1), repeat=m + 1):
                g, o, p = check_profile([((2, 1), e, allocations)])
                gcd_checks += g
                original_gcd_checks += o
                pair_checks += p
                local_profiles += 1
    multi = [
        ((2, 1), 4, (2, 0, 1, 3, 4)),
        ((3, 2), 3, (1, 0, 2, 3, 1)),
        ((5, 2), 2, (0, 0, 1, 2, 1)),
    ]
    g, o, p = check_profile(multi)
    gcd_checks += g
    original_gcd_checks += o
    pair_checks += p
    anchor_checks, residual_checks = check_anchor_invariance(multi)
    moments = check_first_moment()
    hankel_checks, strict = check_hankel_overlap()
    print(
        f"passed {local_profiles} local profiles; "
        f"{gcd_checks} multirow Gaussian gcd checks; "
        f"{original_gcd_checks} original-circle gcd checks; "
        f"{pair_checks} exact pair-numerator identities; "
        f"{anchor_checks} anchor product comparisons "
        f"({residual_checks} nonzero pair residuals); "
        f"{moments} first-moment identities; "
        f"{hankel_checks} overlap Hankel divisibilities "
        f"({strict} strict overlap gains)"
    )


if __name__ == "__main__":
    main()
