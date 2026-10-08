#!/usr/bin/env python3
"""Exact certificate for the obstruction to lifting the four-row fixture.

Uses only Python's standard library. The analytic necessity of the invariant
is proved in quadratic_lift_five_row_moment_obstruction.md.
"""

import json
from fractions import Fraction as Q
from pathlib import Path

from check_four_row_real_polynomial_certificate import main as verify_profile


def invariant(a, b, c):
    return 4 * a * c - 3 * b * b + 4 * a**4


def main():
    verify_profile()
    source = Path(__file__).with_name("four_row_real_polynomial_certificate.json")
    data = json.loads(source.read_text())
    roots = [tuple(map(Q, row)) for row in data["roots"]]
    rho = Q(data["radius"])
    assert len(roots) == 15 and roots[-1] == (0, 1)
    assert rho == Q(1, 10**30)
    assert all(abs(x) + abs(y) + 2 * rho < 2 for x, y in roots)

    a = sum(y for x, y in roots)
    b = sum(2 * x * y for x, y in roots)
    c = sum(3 * x * x * y - y**3 for x, y in roots)
    errors = (14 * rho, 56 * rho, 168 * rho)
    assert all(abs(value) + error < 3
               for value, error in zip((a, b, c), errors))

    e = invariant(a, b, c)
    error = 444 * errors[0] + 18 * errors[1] + 12 * errors[2]
    assert error == 9240 * rho
    assert e - error > 42

    # Finite rational checks supplement, rather than replace, the algebraic
    # identities in the proof. Check the elimination at its forced candidate.
    x, y = b / (2 * a), -a
    assert a + y == 0 and b + 2 * x * y == 0
    assert 4 * a * (c + 3 * x * x * y - y**3) == e
    for scale, shift in ((Q(2), Q(3)), (Q(-3), Q(5, 7)),
                         (Q(7, 11), Q(-13, 17))):
        new_a = scale * a
        new_b = scale**2 * b + 2 * scale * shift * a
        new_c = (scale**3 * c + 3 * scale**2 * shift * b
                 + 3 * scale * shift**2 * a)
        assert invariant(new_a, new_b, new_c) == scale**4 * e
    assert invariant(-a, -b, -c) == e

    # Exact members of the whole vertical three-row family. The formula
    # E=-48*X^2 is proved symbolically in the note; these retain its actual
    # seven-root coordinates and nondegeneracy conditions.
    for u, v in ((Q(1, 2), Q(1, 4)), (Q(1, 3), Q(1, 5))):
        denominator = 1 - u*u - v*v
        s = (1 + u*u + v*v) / denominator
        xx, yy = 2*u / denominator, 2*v / denominator
        parameters = (s*s, -xx*xx, -yy*yy)
        horizontal = s*xx*yy
        assert sum(parameters) == 1
        assert len(set(parameters)) == 3
        assert all(t not in (0, 1, -1) for t in parameters)
        family = [(Q(0), Q(1))]
        family += [(horizontal, t-1) for t in parameters]
        family += [(horizontal*(1+1/t), t) for t in parameters]
        assert len(set(family + [(x, -y) for x, y in family])) == 14
        aa = sum(y for x, y in family)
        bb = sum(2*x*y for x, y in family)
        cc = sum(3*x*x*y-y**3 for x, y in family)
        assert aa == 0 and bb == 4*horizontal
        assert invariant(aa, bb, cc) == -48*horizontal**2 < 0

    print("PASS: the certified four-row profile has E > 42")
    print("No full five-row degree-one extension of any real quadratic pullback")
    print("Center moments:", [float(value) for value in (a, b, c)])
    print(f"Center invariant: {float(e):.15f}")
    print(f"Certified invariant error: {float(error):.3e}")
    print("PASS: vertical three-row family invariant E=-48*X^2")


if __name__ == "__main__":
    main()
