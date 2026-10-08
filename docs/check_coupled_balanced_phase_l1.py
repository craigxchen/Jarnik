#!/usr/bin/env python3
"""Exact small-order checks for the complete balanced-cut phase lattice.

The theorem proved in global documentation gives the all-orders minimum.  This
script checks the formulas at m=6,8,10 without enumerating rational lambda
vectors: integer lattice points are represented by their normalized integer
profiles a_i, and v(S)=sum_{i in S}a_i-sum(a)/2.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, lcm


def cuts(m: int):
    return tuple(combinations(range(m), m // 2))


def values(profile: tuple[int, ...]):
    m = len(profile)
    half = sum(profile) // 2
    return tuple(sum(profile[i] for i in cut) - half for cut in cuts(m))


def lambda_denominator(profile: tuple[int, ...]) -> int:
    """Denominator of the unique sum-zero lambda for this profile."""
    m = len(profile)
    total = sum(profile)
    lambdas = [Fraction(m * a - total, 2 * m) for a in profile]
    return lcm(*(value.denominator for value in lambdas))


def enumerate_short_profiles(m: int, max_level: int = 2):
    """Enumerate normalized profiles with entries in [0,max_level]."""
    best = None
    witnesses = set()
    for profile in product(range(max_level + 1), repeat=m):
        if min(profile) != 0 or sum(profile) % 2:
            continue
        value = values(profile)
        if not any(value):
            continue
        average = Fraction(sum(abs(x) for x in value), comb(m, m // 2))
        if best is None or average < best:
            best = average
            witnesses = {profile}
        elif average == best:
            witnesses.add(profile)
    return best, witnesses


def check_pair_minimizers(m: int) -> None:
    expected = Fraction(m - 2, 2 * (m - 1))
    best, witnesses = enumerate_short_profiles(m)
    assert best == expected
    expected_count = 2 * comb(m, 2)
    assert len(witnesses) == expected_count
    assert {lambda_denominator(p) for p in witnesses} == {m}
    # Every minimizer has exactly two ones or exactly m-2 ones.
    assert all(sum(a == 1 for a in profile) in (2, m - 2)
               for profile in witnesses)
    print(f"m={m}: min average={best}, minimizers={len(witnesses)}, "
          f"lambda denominator(s)={sorted({lambda_denominator(p) for p in witnesses})}")


def check_balanced_cover(m: int) -> None:
    """A two-per-side witness for every balanced cut needs ceil(m/2)+2 rows."""
    target = (m + 1) // 2 + 2
    all_cuts = cuts(m)
    for size in range(4, target):
        for selected in combinations(range(m), size):
            if all(2 <= len(set(selected) & set(cut)) <= size - 2
                   for cut in all_cuts):
                raise AssertionError((m, size, selected))
    selected = tuple(range(target))
    assert all(2 <= len(set(selected) & set(cut)) <= target - 2
               for cut in all_cuts)


if __name__ == "__main__":
    for m in (6, 8, 10):
        check_pair_minimizers(m)
        check_balanced_cover(m)
