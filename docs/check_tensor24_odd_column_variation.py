"""Exact oriented-matroid audit for the degree-22 odd-column fixture.

This checker does not solve the nonlinear odd-moment equations.  It proves
that their five-sign-variation necessary condition is compatible with the
row-difference oriented matroid of the stored code.
"""

import json
from collections import Counter
from itertools import combinations
from math import gcd
from pathlib import Path


ORDER = (6, 17, 8, 14, 4, 9, 18, 21, 0, 7, 12,
         3, 11, 15, 16, 10, 2, 1, 5, 13, 20, 19)
SIGNS = (1, 1, 1, 1, 1, -1, 1, 1, 1, -1, -1,
         -1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1)


def determinant(a):
    """Integer determinant by fraction-free Bareiss elimination."""
    a = [list(row) for row in a]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * value - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
        previous = value
    return sign * a[-1][-1]


def primitive(values):
    divisor = gcd(*(abs(x) for x in values if x))
    out = tuple(x // divisor for x in values)
    if next(x for x in out if x) < 0:
        out = tuple(-x for x in out)
    return out


def variation(c):
    previous = 0
    count = 0
    for j in ORDER:
        value = SIGNS[j] * c[j]
        if not value:
            continue
        value = 1 if value > 0 else -1
        if previous and value != previous:
            count += 1
        previous = value
    return count


def main():
    fixture = json.loads(Path(__file__).with_name(
        "tensor24_odd_column_certificate_gap_fixture.json").read_text())
    s = fixture["retained_matrix"]

    # Scaling by two does not change the row space or any signs.
    b = [[(s[i][j] - s[7][j]) // 2 for j in range(22)]
         for i in range(7)]
    assert sorted(ORDER) == list(range(22))
    assert len(SIGNS) == 22 and set(SIGNS) == {-1, 1}
    pivots = fixture["pivot_columns"]
    assert determinant([[b[i][j] for j in pivots] for i in range(7)])

    # Every row-space cocircuit has a rank-six zero-column flat.  Choosing a
    # basis of that flat gives its normal by the seven signed 6x6 minors.
    cocircuits = set()
    rank_six_subsets = 0
    for zero_basis in combinations(range(22), 6):
        a = [[b[i][j] for i in range(7)] for j in zero_basis]
        normal = []
        for deleted in range(7):
            minor = [row[:deleted] + row[deleted + 1:] for row in a]
            normal.append((-1) ** deleted * determinant(minor))
        if not any(normal):
            continue
        rank_six_subsets += 1
        assert all(sum(row[i] * normal[i] for i in range(7)) == 0
                   for row in a)
        c = [sum(normal[i] * b[i][j] for i in range(7))
             for j in range(22)]
        cocircuits.add(primitive(c))

    support_profile = Counter(sum(x != 0 for x in c) for c in cocircuits)
    variation_profile = Counter(variation(c) for c in cocircuits)
    assert rank_six_subsets == 69440
    assert len(cocircuits) == 10155
    assert support_profile == {
        7: 1, 10: 1, 11: 44, 12: 105, 13: 154,
        14: 639, 15: 2086, 16: 7125,
    }
    assert min(variation_profile) == 5

    support_seven = next(c for c in cocircuits
                         if sum(x != 0 for x in c) == 7)
    assert support_seven == (
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, 0, 0, 0, -2,
        -3, -1, 0, -2, 0, -2,
    )
    print("PASS: 10,155 exact row-space cocircuits;")
    print("minimum support 7; compatible minimum sign variation 5.")


if __name__ == "__main__":
    main()
