#!/usr/bin/env python3
"""Exhaust the odd-subset row selections in two four-order cyclotomic families.

This dependency-free search enumerates all 255 nonempty subsets of the eight
odd-cardinality masks of four orders. It computes the exact divisor-union
totient weights and actual orientation width. The search is exhaustive only
for these eight candidate rows at each displayed order fixture; it is not a
classification of other cyclotomic supports.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd, isqrt


ODD_MASKS = tuple(m for m in range(16) if bin(m).count("1") % 2 == 1)


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        if n % p == 0:
            a = 0
            while n % p == 0:
                n //= p
                a += 1
            out[p] = a
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = 1
    return out


def divisor_phi_table(n: int) -> dict[int, int]:
    """Return {d: phi(d)} for all d|n, by prime-power multiplication."""
    table = {1: 1}
    for p, exponent in factor(n).items():
        old = tuple(table.items())
        powers: list[tuple[int, int]] = []
        pe = phipe = 1
        for a in range(1, exponent + 1):
            pe *= p
            phipe = p - 1 if a == 1 else phipe * p
            powers.append((pe, phipe))
        for d, phid in old:
            for pe, phipe in powers:
                table[d * pe] = phid * phipe
    return table


def order_family(k: int) -> tuple[int, int, int, int]:
    assert k >= 19 and k % 2 == 1 and k % 3 in (0, 1)
    x, y, z, w = k - 2, k, k + 2, k + 4
    assert all(gcd(a, b) == 1 for a, b in combinations((x, y, z, w), 2))
    return x * z * w, y * y * w, y * z * z, y * z * w


def squarefree_fixture() -> tuple[int, int, int, int]:
    # The p-close construction from cyclotomic_nullity_four_counterexample.md.
    p1, p2, p3 = 181, 193, 199
    factors = (179, 191, 197)
    assert all(
        all(q % d for d in range(2, isqrt(q) + 1))
        for q in (p1, p2, p3, *factors)
    )
    q = p1 * p2 * p3
    ns = (179 * q // p1, 191 * q // p2, 197 * q // p3, q)
    assert ns == (6874853, 6879629, 6881801, 6951667)
    return ns


def width_summary(ns: tuple[int, ...]) -> dict[str, object]:
    h = min(ns)
    assert len(set(ns)) == 4 and max(ns) < 3 * h
    tables = [divisor_phi_table(n) for n in ns]
    phi_by_e: dict[int, int] = {}
    for table in tables:
        for e, value in table.items():
            if e in phi_by_e:
                assert phi_by_e[e] == value
            phi_by_e[e] = value
    es = tuple(sorted(phi_by_e))
    assert tuple(e for e in es if e >= h) == tuple(sorted(ns))

    # A row mask selects one of the four indicator vectors 1_{e|n_j}; its
    # coordinate is the number of selected orders divisible by e.
    values = {
        mask: tuple(sum((mask >> j) & 1 for j, n in enumerate(ns) if n % e == 0)
                    for e in es)
        for mask in ODD_MASKS
    }
    high_positions = tuple(es.index(n) for n in ns)
    for mask in ODD_MASKS:
        assert tuple(values[mask][pos] for pos in high_positions) == tuple(
            (mask >> j) & 1 for j in range(4)
        )
    by_size: dict[int, dict[str, object]] = {}
    for size in range(1, len(ODD_MASKS) + 1):
        eligible: list[tuple[int, tuple[int, ...]]] = []
        for selection in combinations(ODD_MASKS, size):
            # Every high coordinate must vary. This makes n_1=h appear in a
            # pair difference, so the minimum contact is exactly h; all low
            # Möbius coordinates vanish for every pair.
            if not all(len({values[m][pos] for m in selection}) > 1
                       for pos in high_positions):
                continue
            widths = tuple(max(values[m][j] for m in selection) -
                           min(values[m][j] for m in selection)
                           for j in range(len(es)))
            base_width = widths[es.index(1)]
            assert base_width % 2 == 0
            L = sum(phi_by_e[e] * widths[j] for j, e in enumerate(es))
            eligible.append((L, selection))
        if eligible:
            min_L = min(L for L, _ in eligible)
            min_sets = tuple(selection for L, selection in eligible if L == min_L)
            by_size[size] = {
                "eligible": len(eligible),
                "min_L": min_L,
                "min_sets": min_sets,
                "at_most_4h": sum(L <= 4 * h for L, _ in eligible),
                "strictly_below_4h": sum(L < 4 * h for L, _ in eligible),
            }
        else:
            by_size[size] = {
                "eligible": 0,
                "min_L": None,
                "min_sets": (),
                "at_most_4h": 0,
                "strictly_below_4h": 0,
            }

    return {
        "h": h,
        "orders": ns,
        "divisor_count": len(es),
        "phi_union": sum(phi_by_e.values()),
        "four_h": 4 * h,
        "by_size": by_size,
    }


def check_k19() -> None:
    summary = width_summary(order_family(19))
    assert summary["orders"] == (8211, 8303, 8379, 9177)
    assert summary["divisor_count"] == 36
    assert summary["phi_union"] == 32751
    assert summary["h"] == 8211 and summary["four_h"] == 32844
    expected = {
        5: (56, 32812, 2),
        6: (28, 33586, 0),
        7: (8, 34046, 0),
        8: (1, 34068, 0),
    }
    for size, (count, minimum, fitting) in expected.items():
        result = summary["by_size"][size]
        assert result["eligible"] == count
        assert result["min_L"] == minimum
        assert result["at_most_4h"] == fitting
        assert result["strictly_below_4h"] == fitting
    assert summary["by_size"][5]["min_sets"] == (
        (1, 2, 4, 7, 8), (7, 8, 11, 13, 14)
    )
    print("k=19:", summary)


def check_squarefree_fixture() -> None:
    summary = width_summary(squarefree_fixture())
    assert summary["divisor_count"] == 20
    assert summary["phi_union"] == 27478591
    assert summary["h"] == 6874853 and summary["four_h"] == 27499412
    expected_minima = {
        5: 27479162,
        6: 27549542,
        7: 27587750,
        8: 27587948,
    }
    for size, minimum in expected_minima.items():
        result = summary["by_size"][size]
        assert result["min_L"] == minimum
        assert result["at_most_4h"] == (2 if size == 5 else 0)
        assert result["strictly_below_4h"] == (2 if size == 5 else 0)
    print("squarefree fixture:", summary)


def check_samples() -> None:
    # Finite samples in the same coprime odd-k family, not an assertion for
    # all k. Each listed k meets the hypotheses of order_family.
    expected_size6_gap = {
        19: -742,
        21: -858,
        31: -1558,
        51: -3558,
        61: -4858,
        81: -8058,
        91: -9958,
        201: -44058,
    }
    for k, gap in expected_size6_gap.items():
        summary = width_summary(order_family(k))
        actual = summary["four_h"] - summary["by_size"][6]["min_L"]
        assert actual == gap
        print(f"sample k={k}: best size-6 margin 4h-L={actual}")


def main() -> None:
    assert len(ODD_MASKS) == 8
    check_k19()
    check_squarefree_fixture()
    check_samples()
    print("All exact odd-subset search checks passed.")


if __name__ == "__main__":
    main()
