"""The two excess interpolation equations on an eight-point affine plane."""
import sympy as sp

from ii_min_class import coordinates

a, b = sp.symbols('a b')


def interpolation(ids):
    h = coordinates(ids)
    assert h.shape == (8, 3) and [0, 0, 1] in h.tolist()
    base = sp.Matrix([a, b, 1])
    heights = h*base
    design = sp.Matrix([[x*x, y*y, z*z, x*y, x*z, y*z]
                        for x, y, z in h.tolist()])
    pivots = list(design.T.rref()[1])
    assert len(pivots) == 6
    coeffs = design[pivots, :].inv() * sp.Matrix([heights[i]**4 for i in pivots])
    equations = [sp.factor((design.row(i)*coeffs)[0] - heights[i]**4)
                 for i in range(h.rows) if i not in pivots]
    return h, heights, equations
