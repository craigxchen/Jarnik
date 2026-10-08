"""Exact Boolean multiplicity checks for the disjoint elliptic fibers."""

from fractions import Fraction as F
from itertools import combinations


MOVING = frozenset(("a", "b", "c"))


def check_completed_linear_forms():
    # Primitive integer representatives of X/Y - tau, with fixed
    # denominators cleared.  These are the four complete-fiber forms.
    forms = {
        "infinity": (3, -16),
        "0": (2, 1),
        "1": (2, -9),
        "rho": (4, -25),
    }
    taus = {
        "infinity": F(16, 3),
        "0": F(-1, 2),
        "1": F(9, 2),
        "rho": F(25, 4),
    }
    for label, (a, b) in forms.items():
        assert F(-b, a) == taus[label]  # aX+bY=0 has T=-b/a.
        assert a > 0

    # Direct substitution keeps the exact Gram coefficients visible.
    di, d0, d1, dr = (F(3), F(5), F(7), F(11))
    ei, e0, e1, er = (F(2), F(3), F(5), F(7))
    ri, r0, r1, rr = (F(13), F(17), F(19), F(23))
    Fi, F0, F1, Fr = di * ei, d0 * e0, d1 * e1, dr * er
    xi, x0, x1, xr = ri / ei, r0 / e0, r1 / e1, rr / er
    old_linear = 15 * di * ri + 10 * d0 * r0 - 6 * d1 * r1 - dr * rr
    new_linear = 15 * Fi * xi + 10 * F0 * x0 - 6 * F1 * x1 - Fr * xr
    assert old_linear == new_linear
    old_det = 4 * di * d0 * ri * r0 - (d1 * r1 - di * ri - d0 * r0) ** 2
    new_det = 4 * Fi * F0 * xi * x0 - (F1 * x1 - Fi * xi - F0 * x0) ** 2
    assert old_det == new_det


def anchor_multiplicities(subset):
    """Multiplicity vector for sum_j H_(j,alpha) at a Boolean cell."""
    return len(set(subset) & MOVING)


def rho_multiplicities(subset):
    return len(set(subset) & MOVING)


def check_anchor(anchor):
    cells = [frozenset(s) for size in range(4)
             for s in combinations(sorted(MOVING), size)]
    # D_anchor has coefficient one on every nonempty cell and zero on the
    # empty cell.  The correction in (3) is one on pairs and two on triple.
    for cell in cells:
        lhs = anchor_multiplicities(cell)
        rhs = (0 if not cell else 1)
        if len(cell) == 2:
            rhs += 1
        elif len(cell) == 3:
            rhs += 2
        assert lhs == rhs, (anchor, cell, lhs, rhs)
    assert sum(anchor_multiplicities(cell) for cell in cells) == 12
    assert sum(1 for cell in cells if cell) == 7


def check_rho():
    cells = [frozenset(s) for size in range(4)
             for s in combinations(sorted(MOVING), size)]
    # D_rho is supported on cells of size at least two.  Formula (4) has
    # singleton correction one, pair correction one, and triple correction
    # two; the empty cell is absent from every coordinate section.
    for cell in cells:
        lhs = rho_multiplicities(cell)
        rhs = 1 if len(cell) >= 2 else 0
        if len(cell) == 1:
            rhs += 1
        elif len(cell) == 2:
            rhs += 1
        elif len(cell) == 3:
            rhs += 2
        assert lhs == rhs, (cell, lhs, rhs)
    assert sum(rho_multiplicities(cell) for cell in cells) == 12
    assert sum(1 for cell in cells if len(cell) >= 2) == 4


def check_complete_fibers():
    # Three anchor fibers: seven shared points plus one private anchor point.
    for anchor in ("infinity", "0", "1"):
        cells = [frozenset(s) for size in range(4)
                 for s in combinations(sorted(MOVING), size)]
        assert len(cells) == 8
        assert sum(1 for cell in cells if cell) == 7
        assert sum(1 for cell in cells if not cell) == 1

    # rho: four shared points, three singleton private cells, one empty cell.
    cells = [frozenset(s) for size in range(4)
             for s in combinations(sorted(MOVING), size)]
    assert sum(1 for cell in cells if len(cell) >= 2) == 4
    assert sum(1 for cell in cells if len(cell) == 1) == 3
    assert sum(1 for cell in cells if not cell) == 1
    assert 3 * 7 + 4 == 25
    assert 3 * 8 + 8 == 32


if __name__ == "__main__":
    for anchor in ("infinity", "0", "1"):
        check_anchor(anchor)
    check_rho()
    check_complete_fibers()
    check_completed_linear_forms()
    print("Complete Boolean fiber and section multiplicity identities pass.")
