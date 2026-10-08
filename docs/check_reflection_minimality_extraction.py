#!/usr/bin/env python3
"""Exact checks for reflection minimality under row restriction."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations


def width(values: tuple[int, ...], rows: tuple[int, ...]) -> int:
    selected = [values[r] for r in rows]
    return max(selected) - min(selected)


def delta(values: tuple[int, ...], rows: tuple[int, ...],
          a: int, b: int, ell: int) -> int:
    old = width(values, rows)
    new_value = values[a] + values[b] - values[ell]
    changed = [values[r] for r in rows if r != ell] + [new_value]
    return max(changed) - min(changed) - old


def reflection_triples(rows: tuple[int, ...]):
    for ell in rows:
        retained = [r for r in rows if r != ell]
        for pos, a in enumerate(retained):
            for b in retained[pos:]:
                yield a, b, ell


def ratio(primes: tuple[int, ...], changes: tuple[int, ...]) -> Fraction:
    answer = Fraction(1)
    for prime, change in zip(primes, changes):
        answer *= Fraction(prime ** max(change, 0), prime ** max(-change, 0))
    return answer


def uniform_gcd_counts(M: int, k: int) -> None:
    rows = tuple(range(M))
    subset = set(range(k))
    ell = 0
    retained = subset - {ell}
    outside = set(rows) - subset

    # Choose the side containing the indicated retained set. Count extra
    # cuts that become part of gcd(rows in I\{ell}) after restriction.
    extra_retained_private = 0
    subset_common = 0
    full = (1 << M) - 1
    cuts = {min(mask, full ^ mask) for mask in range(1, full)}
    for key in cuts:
        sides = ({r for r in rows if (key >> r) & 1},
                 {r for r in rows if not ((key >> r) & 1)})
        if any(retained <= side and (other & outside)
               for side, other in (sides, sides[::-1])):
            extra_retained_private += 1
        if any(subset <= side and bool(other & outside)
               for side, other in (sides, sides[::-1])):
            subset_common += 1

    assert extra_retained_private == 2 ** (M-k+1) - 2
    assert subset_common == 2 ** (M-k) - 1
    raw_retained_gcd = 1 + extra_retained_private
    intrinsic_retained_gcd = raw_retained_gcd - subset_common
    assert raw_retained_gcd == 2 ** (M-k+1) - 1
    assert intrinsic_retained_gcd == 2 ** (M-k)


def local_minimality_counterexample() -> None:
    # Columns are valuations at Gaussian primes above 5, 13, and 17.
    columns = (
        (0, 2, 2, 0, 1),
        (2, 2, 0, 0, 2),
        (1, 1, 0, 2, 0),
    )
    primes = (5, 13, 17)
    full_rows = tuple(range(5))

    # Every one-coordinate reflection of the full five-row tuple preserves
    # or increases its primitive common norm. Repeated anchors a=b are
    # included.
    for a, b, ell in reflection_triples(full_rows):
        changes = tuple(delta(column, full_rows, a, b, ell)
                        for column in columns)
        assert ratio(primes, changes) >= 1

    # After extracting rows 0,1,3,4, the reflection replacing row 3 by
    # z_0 z_1/z_3 decreases the primitive common norm by exactly 17.
    subset_rows = (0, 1, 3, 4)
    changes = tuple(delta(column, subset_rows, 0, 1, 3)
                    for column in columns)
    assert changes == (0, 0, -1)
    assert ratio(primes, changes) == Fraction(1, 17)


def main() -> None:
    for M in range(4, 9):
        for k in range(3, M):
            uniform_gcd_counts(M, k)
    local_minimality_counterexample()
    print("PASS: uniform-cut gcd amplification holds for 3<=k<M<=8.")
    print("PASS: the five-row tuple is reflection-minimal, but its four-row restriction drops by 17.")


if __name__ == "__main__":
    main()
