#!/usr/bin/env python3
"""Finite exact check for affine_full_fair_to_half_symmetric_obstruction.md."""

from itertools import combinations, permutations, product


def extrema_have_two_each(values):
    lo = min(values)
    hi = max(values)
    return lo == hi or (values.count(lo) >= 2 and values.count(hi) >= 2)


def partition_of_nonconstant(values):
    if len(set(values)) == 1:
        return None
    lo = min(values)
    hi = max(values)
    if values.count(lo) != 2 or values.count(hi) != 2:
        return False
    low = frozenset(i for i, value in enumerate(values) if value == lo)
    high = frozenset(set(range(4)) - low)
    return frozenset((low, high))


def main():
    # All bounded integer column vectors that can occur after the singleton
    # test: constants and permutations of (u,u,v,v).  The row-sum-one
    # condition is imposed across source columns, not down a column.
    vectors = set()
    for value in range(-2, 3):
        vectors.add((value, value, value, value))
    for u in range(-2, 3):
        for v in range(-2, 3):
            if u != v:
                vectors.update(permutations((u, u, v, v)))
    vectors = sorted(vectors)
    nonconstant = []
    for v in vectors:
        part = partition_of_nonconstant(v)
        if part not in (None, False):
            nonconstant.append((v, part))

    # Every pair of nonconstant vectors from different partitions has a
    # forbidden singleton extremum after summation.
    checked = 0
    for (u, part_u), (v, part_v) in combinations(nonconstant, 2):
        if part_u == part_v:
            continue
        summed = tuple(a + b for a, b in zip(u, v))
        assert not extrema_have_two_each(summed)
        checked += 1
    assert checked > 0

    # Exhaustively check the implication for bounded collections of columns:
    # if every singleton and every pair sum is admissible, all nonconstant
    # columns share one partition.
    for cols in combinations(vectors, 4):
        if not all(extrema_have_two_each(v) for v in cols):
            continue
        if not all(
            extrema_have_two_each(tuple(a + b for a, b in zip(u, v)))
            for u, v in combinations(cols, 2)
        ):
            continue
        parts = {
            partition_of_nonconstant(v)
            for v in cols
            if partition_of_nonconstant(v) not in (None, False)
        }
        assert len(parts) <= 1

    # Exact radius-range bookkeeping on one sample exponent column.
    sample = ((-1, 0, 1, 1), (0, -1, 1, 1), (0, 0, 0, 1), (0, 0, 1, 0))
    exponents = tuple(sum(row[j] for row in sample) for j in range(4))
    assert exponents == (-1, -1, 3, 3)
    assert max(exponents) - min(exponents) == 4
    residual = tuple(x - min(exponents) for x in exponents)
    residual_bar = tuple(max(exponents) - x for x in exponents)
    assert residual == (0, 0, 4, 4)
    assert residual_bar == (4, 4, 0, 0)

    print("affine full-fair half-symmetric obstruction: PASS")
    print("checked cross-partition pairs:", checked)


if __name__ == "__main__":
    main()
