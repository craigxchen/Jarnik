#!/usr/bin/env python3
"""Exact bounded checks for the whole-clique affine translation identity."""

import sys
from itertools import combinations
from math import lcm

sys.path.insert(0, "docs")
from check_integer_cotangent_normalization import check_clique, edge_quotient


def shift_data(xs, L, t):
    ys = tuple(x + t for x in xs)
    for i, j in combinations(range(len(xs)), 2):
        d = xs[i] - xs[j]
        c = edge_quotient(xs[i], xs[j], L)
        assert edge_quotient(ys[i], ys[j], L) == c + t * (xs[i] + xs[j] + t) // d
    return ys


def canonical_shift(xs):
    return lcm(*(abs(a - b) for a, b in combinations(xs, 2)))


if __name__ == "__main__":
    cases = [((8, 9, 10, 12), 6), ((-18, -6, -3, 0, 2, 6, 12), 6)]
    for xs, L in cases:
        N = check_clique(xs, L)
        D = canonical_shift(xs)
        print(f"X={xs}, L={L}, D={D}, N={N}")
        for t in (D, -D):
            ys = shift_data(xs, L, t)
            N1 = check_clique(ys, L)
            print(f"  t={t}: N'={N1}, ratio={N1}/{N}, positive={min(ys)>0}")
    print("PASS: exact whole-clique translation closure and lcm updates")
