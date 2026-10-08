#!/usr/bin/env python3
"""Exact tests for bounded extraction of the global quartet content.

The rows are equal-norm Gaussian integers.  For each label ``i`` we use one
source block kappa_i on the majority rows and its conjugate on row i.  Thus
every source prime has a singleton residue class, with a different singleton
for each block.  The script compares the global quartet gcd of all rows with
the gcd after retaining a subset of rows.

This is a finite obstruction to a bound independent of m.  It is a valuation
profile obstruction and makes no claim that the resulting rows fit an
endpoint arc of length C*sqrt(R).
"""

from __future__ import annotations

from itertools import combinations
from math import comb, gcd, lcm, prod

from check_quartet_matching_gcd_cut_budget import (ggcd, matching_content,
                                                    mul, norm, sub)
from check_global_matching_triangle_content_normalization import (
    factors, gaussian_prime, vg,
)


Gaussian = tuple[int, int]


def conjugate(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def source_blocks(m: int, q: int = 1) -> tuple[Gaussian, ...]:
    M = lcm(*(j * j - k * k for j in range(1, m + 1)
              for k in range(1, j))) if m >= 2 else 1
    return tuple((2 * j * M * q, 1) for j in range(1, m + 1))


def singleton_rows(m: int, q: int = 1) -> tuple[tuple[Gaussian, ...], tuple[Gaussian, ...]]:
    """Return rows and source blocks for the singleton profile."""
    blocks = source_blocks(m, q)
    block_norms = tuple(norm(block) for block in blocks)
    assert all(gcd(a, b) == 1 for a, b in combinations(block_norms, 2))
    rows: list[Gaussian] = []
    for outlier in range(m):
        row = (1, 0)
        for j, block in enumerate(blocks):
            row = mul(row, conjugate(block) if j == outlier else block)
        rows.append(row)
    assert len(set(rows)) == m
    assert len({norm(row) for row in rows}) == 1
    row_gcd = (0, 0)
    for row in rows:
        row_gcd = ggcd(row_gcd, row)
    assert norm(row_gcd) == 1
    return tuple(rows), blocks


def global_content(rows: tuple[Gaussian, ...]) -> Gaussian:
    """Gaussian gcd of all two-edge matching products on quartets."""
    answer = (0, 0)
    for inds in combinations(range(len(rows)), 4):
        answer = ggcd(answer, matching_content(tuple(rows[i] for i in inds)))
    return answer


def preserves(full: Gaussian, candidate: Gaussian) -> bool:
    """The candidate gcd is the same ideal as full (units are immaterial)."""
    # A subset has fewer quartet factors, hence its gcd is divisible by full.
    # Equality of norms is therefore equivalent to equality of ideals.
    return norm(candidate) == norm(full)


def minimum_witness(rows: tuple[Gaussian, ...]) -> tuple[int, tuple[int, ...]]:
    full = global_content(rows)
    for size in range(4, len(rows) + 1):
        for inds in combinations(range(len(rows)), size):
            candidate = global_content(tuple(rows[i] for i in inds))
            if preserves(full, candidate):
                return size, inds
    raise AssertionError("the full row set must preserve its own content")


def audit_singleton_family(m: int) -> None:
    rows, blocks = singleton_rows(m)
    block_norms = tuple(norm(block) for block in blocks)
    full = global_content(rows)
    size, witness = minimum_witness(rows)
    # Each block is a source factor in the full G, and its unique outlier must
    # be retained.  The exact enumeration below certifies the conclusion.
    assert size == m and witness == tuple(range(m))
    assert norm(full) % prod(block_norms) == 0
    # For every rational source factor p of block j, the selected Gaussian
    # orientation has exponent e on the full G.  Omitting j makes every
    # selected difference divisible by that orientation, so its quartet gcd
    # has exponent at least 2e.  The finite valuation check records this
    # strict jump directly.
    for j, block_norm in enumerate(block_norms):
        majority = [i for i in range(m) if i != j]
        for p in factors(block_norm):
            pi = gaussian_prime(p)
            if vg(blocks[j], pi) == 0:
                pi = (pi[0], -pi[1])
            exponent = vg(blocks[j], pi)
            assert exponent > 0 and vg(full, pi) == exponent
            assert min(vg(sub(rows[b], rows[a]), pi)
                       for a, b in combinations(majority, 2)) == exponent
    print(f"m={m}: NormG={norm(full)}, block-product={prod(block_norms)}, "
          f"minimum witness size={size}, witness={witness}")


def audit_cut_requirement() -> None:
    """Check the exact singleton requirement in a symbolic cut model.

    A quartet with the singleton outlier and three majority rows has one
    cross edge and one majority edge, so its source valuation is e+r.  A
    quartet omitting the outlier has two majority edges and cannot attain the
    e+r valuation when e>r.  This is the only combinatorial fact used by the
    family above; the exact Gaussian rows provide the finite certificate.
    """
    for m in range(4, 13):
        for outlier in range(m):
            majority = [i for i in range(m) if i != outlier]
            assert len(majority) >= 3
            assert len(list(combinations(majority, 3))) == comb(m - 1, 3)


def minimum_balanced_cover(m: int) -> tuple[int, tuple[int, ...]]:
    """Smallest T meeting every m/2-sized cut in two labels each side."""
    n = m // 2
    cuts = tuple(set(cut) for cut in combinations(range(m), n))
    for size in range(4, m + 1):
        for selected in combinations(range(m), size):
            if all(2 <= len(set(selected) & cut) <= size - 2 for cut in cuts):
                return size, selected
    raise AssertionError("the full row set covers every balanced cut")


def audit_balanced_cover() -> None:
    for m in range(4, 11):
        size, selected = minimum_balanced_cover(m)
        expected = m // 2 + 2 + (m % 2)
        assert size == expected
        assert selected == tuple(range(size))


def audit_compact_fixture() -> None:
    blocks = tuple((a, 1) for a in (4, 6, 10, 14, 16))
    rows = []
    for outlier in range(5):
        row = (1, 0)
        for j, block in enumerate(blocks):
            row = mul(row, conjugate(block) if j == outlier else block)
        rows.append(row)
    rows = tuple(rows)
    assert rows == (
        (56046, 8675), (53954, 17475), (51210, 24371),
        (49746, 27235), (49254, 28115),
    )
    assert norm(global_content(rows)) == 823400893696
    assert minimum_witness(rows) == (5, (0, 1, 2, 3, 4))


if __name__ == "__main__":
    audit_cut_requirement()
    audit_balanced_cover()
    audit_compact_fixture()
    for m in range(4, 9):
        audit_singleton_family(m)
