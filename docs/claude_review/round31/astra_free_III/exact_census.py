"""Enumerate pair-admissible six-row configurations for D=8, c=0.

Build unordered sets of rows modulo row permutation, coordinate permutation
and complementing all rows. Apply the saturated mass test after enumeration.
"""
from functools import lru_cache
from itertools import permutations, product

import numpy as np
import sympy as sp


def byte_rows(columns, n):
    return [sum(((col >> i) & 1) << j for j, col in enumerate(columns)) for i in range(n)]


def row_matrix(columns, n):
    return sp.Matrix([[(row >> j) & 1 for j in range(8)] for row in byte_rows(columns, n)])


@lru_cache(None)
def tables(n):
    return np.array([
        [sum(((col >> i) & 1) << j for j, i in enumerate(perm)) for col in range(2**n)]
        for perm in permutations(range(n))
    ], dtype=np.int64)


def canonical(cols, n):
    # Encode eight sorted column masks as base-2**n digits. Every
    # code has at most 48 bits, so the integer dot products are exact.
    transformed = np.sort(tables(n)[:, cols], axis=1)
    weights = 1 << np.arange(0, 8*n, n, dtype=np.int64)
    flipped = (2**n - 1) - transformed[:, ::-1]
    code = int(min(np.min(transformed @ weights), np.min(flipped @ weights)))
    return tuple((code >> (i*n)) % (1 << n) for i in range(8))


def pair_ok(x, y):
    return (x & ~y).bit_count() >= 2 and (y & ~x).bit_count() >= 2


def cliques(last=6):
    pool = {(0,) * 8}
    for n in range(1, last+1):
        new = set()
        seen = set()
        for cols in pool:
            rows = byte_rows(cols, n-1)
            for row in range(256):
                if not all(pair_ok(r, row) for r in rows):
                    continue
                key = tuple(sorted(col + (((row >> j) & 1) << (n-1))
                                   for j, col in enumerate(cols)))
                if key in seen:
                    continue
                seen.add(key)
                new.add(canonical(key, n))
        print(f'{n} rows: {len(new)} pair-admissible classes')
        pool = new
    return pool


def saturated_filter(classes):
    # Negative-sum tests are represented by sign reversal. This finite
    # box gives a necessary condition even if other forbidden vectors exist.
    forbidden = [v for v in product(range(-2, 3), repeat=8)
                 if any(v) and sum(max(z, 0) for z in v) <= 1]
    tests = np.array(forbidden, dtype=np.int64)
    survivors = []
    rank_five = 0
    for cols in sorted(classes):
        rows = row_matrix(cols, 6)
        differences = rows[1:, :] - sp.ones(5, 1) * rows[0, :]
        kernel = differences.nullspace()
        if len(kernel) != 3:
            continue
        rank_five += 1
        basis = sp.Matrix.hstack(*kernel)
        integral = basis * sp.ilcm(*[v.q for v in basis])
        assert 16 * max(abs(int(v)) for v in integral) < 2**63
        if np.any(np.all(tests @ np.array(integral, dtype=np.int64) == 0, axis=1)):
            continue
        survivors.append(cols)
    print(f'rank five: {rank_five}; surviving test: {len(survivors)}')
    return survivors


if __name__ == '__main__':
    saturated_filter(cliques())
