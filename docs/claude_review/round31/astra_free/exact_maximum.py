"""Enumerate the D=8 extensions in an affine four-plane of five given rows.

Run python3 exact_maximum.py. This supplements the rank-five exclusion in
exact_plane.py. See the report for the elementary rank-at-most-three case.
"""
from itertools import combinations

import numpy as np
import sympy as sp

from exact_census import byte_rows, cliques, pair_ok, row_matrix

cube = np.array([[(r >> j) & 1 for j in range(8)] for r in range(256)], dtype=np.int64)


def maximum(codes):
    best = 0
    extremizers = set()
    for cols in codes:
        basis = byte_rows(cols, 5)
        rows = row_matrix(cols, 5)
        orthogonal = (rows[1:, :] - sp.ones(4, 1)*rows[0, :]).nullspace()
        if len(orthogonal) != 4:
            continue
        H = sp.Matrix.hstack(*orthogonal)
        H *= sp.ilcm(*[v.q for v in H])
        assert 8 * max(abs(h) for h in H) < 2**63
        offsets = (cube - cube[basis[0]]) @ np.array(H, dtype=np.int64)
        affine_plane = np.flatnonzero(np.all(offsets == 0, axis=1))
        pool = [int(r) for r in affine_plane if all(pair_ok(int(r), s) for s in basis)]
        for k in range(len(pool)+1):
            for extra in combinations(pool, k):
                if not all(pair_ok(r, s) for r, s in combinations(extra, 2)):
                    continue
                if 5+k > best:
                    best = 5+k
                    extremizers.clear()
                if 5+k == best:
                    extremizers.add(tuple(sorted([*basis, *extra])))
    return best, sorted(extremizers)


if __name__ == '__main__':
    best, sets = maximum(cliques(5))
    assert best == 8 and len(sets) == 3
    for rows in sets:
        print('rank-four maximizer from canonical five-row base', rows)

    perm = [0, 4, 1, 5, 2, 6, 3, 7]
    permuted = tuple(sorted(sum(((r >> j) & 1) << k for j, k in enumerate(perm))
                            for r in sets[0]))
    assert permuted == sets[1]

    rows = sp.Matrix(cube[list(sets[2]), :])
    H = sp.Matrix.hstack(*(rows[1:, :] - sp.ones(7, 1)*rows[0, :]).nullspace())
    assert H[0, :] == H[7, :]
    print('First two classes are equivalent; third requires roots 0 and 7 to coincide.')
