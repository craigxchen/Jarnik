"""Exact nonlinear certificates for the round-30 D=8, c=0 problem.

Run python3 exact_plane.py from this directory. The census is recomputed,
and its classes are matched to the recorded zero-based order for reporting.
Requires SymPy and NumPy. No access to the original repository is needed.
"""
from itertools import combinations

import sympy as sp

from exact_census import canonical, cliques, saturated_filter

# Each number encodes one row, with bit j set when its j-th sign is +1.
ARCHIVED_ROWS = [
    (7, 57, 201, 90, 178, 156), (7, 57, 201, 90, 178, 108),
    (7, 57, 201, 90, 178, 172), (7, 57, 201, 90, 178, 228),
    (7, 57, 201, 90, 162, 116), (7, 57, 201, 90, 156, 240),
    (7, 57, 193, 74, 146, 100), (7, 25, 97, 42, 146, 84),
    (7, 25, 97, 42, 146, 196), (15, 51, 195, 85, 169, 150),
    (15, 51, 195, 85, 169, 156), (15, 51, 195, 85, 150, 108),
    (15, 51, 195, 85, 166, 156), (15, 51, 195, 85, 166, 120),
    (15, 51, 195, 85, 156, 228), (15, 51, 195, 85, 156, 232),
    (15, 51, 85, 169, 150, 202), (15, 51, 85, 169, 150, 228),
    (15, 51, 85, 169, 198, 180), (15, 51, 85, 169, 198, 120),
    (15, 51, 85, 169, 198, 216),
]
a, b = sp.symbols('a b')


def archive_key(rows):
    cols = tuple(sum(((row >> j) & 1) << i for i, row in enumerate(rows))
                 for j in range(8))
    return canonical(cols, 6)


def normalized_basis(rows):
    signs = sp.Matrix([[(row >> j) & 1 for j in range(8)] for row in rows])
    differences = signs[1:, :] - sp.ones(5, 1) * signs[0, :]
    H = sp.Matrix.hstack(*differences.nullspace())
    assert H.shape == (8, 3)
    # The root at this index has y=1 after common positive dilation.
    assert [0, 0, 1] in H.tolist()
    return H


def geometric_certificate(H):
    for inds in combinations(range(8), 2):
        if H[list(inds), :].rank() < 2:
            return 'proportional rows', inds
    lifted = H.row_join(sp.ones(8, 1))
    for inds in combinations(range(8), 4):
        if lifted[list(inds), :].rank() <= 2:
            return 'four collinear rows', inds
    return None


def plane_certificate(H):
    u = sp.Matrix([a, b, 1])
    heights = H * u
    h = [H.row(j).T for j in range(8)]
    # Coefficients in the basis of six quadratic monomials on Q^3.
    design = sp.Matrix([[x*x, y*y, z*z, x*y, x*z, y*z] for x, y, z in h])
    pivots = list(design.T.rref()[1])
    assert len(pivots) == 6
    coefficients = design[pivots, :].inv() * sp.Matrix([heights[j]**4 for j in pivots])
    xx, yy, zz, xy, xz, yz = coefficients
    Q = sp.Matrix([[xx, xy / 2, xz / 2],
                   [xy / 2, yy, yz / 2],
                   [xz / 2, yz / 2, zz]])
    equations = [(h[j].T*Q*h[j])[0] - heights[j]**4 for j in range(8) if j not in pivots]
    # Restrict Q to u-perp, where Q equals 3 vv^T.
    F = sp.Matrix([[1, 0, -a], [0, 1, -b]])
    restriction = (F*Q*F.T).applyfunc(sp.expand)
    equations.append(restriction.det())
    G = sp.groebner(equations, b, a, order='lex', domain=sp.QQ)

    def reduce(p):
        return G.reduce(sp.expand(p))[1]

    if reduce(1) == 0:
        return 'unit ideal', ()
    # At every admissible solution the product is positive. Reducing after
    # each multiplication preserves its residue class and keeps degrees small.
    product = sp.Integer(1)
    for value in heights:
        product = reduce(product * value)
        if product == 0:
            return 'vanishing height product', ()
    for j, k in combinations(range(8), 2):
        t = h[j] * heights[k] - h[k] * heights[j]
        squared_distance = ((heights[j]-heights[k]) * heights[j] * heights[k])**2 + (t.T*Q*t)[0]/3
        product = reduce(product * reduce(squared_distance))
        if product == 0:
            return 'vanishing distance product', (j, k)
    raise ArithmeticError('No certificate found')


if __name__ == '__main__':
    survivors = saturated_filter(cliques())
    archive = {archive_key(rows) for rows in ARCHIVED_ROWS}
    assert len(archive) == len(ARCHIVED_ROWS) == 21 and archive == set(survivors)
    print('Independent census matches the 21 archived classes.')
    for index, rows in enumerate(ARCHIVED_ROWS):
        H = normalized_basis(rows)
        kind, indices = geometric_certificate(H) or plane_certificate(H)
        print(index, kind, indices)
