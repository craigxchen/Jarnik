"""Exact local and torsion certificates for the two-root obstruction."""

from fractions import Fraction as F
from itertools import product
from math import gcd


def square_residues(modulus):
    return {x * x % modulus for x in range(modulus)}


def check_local_certificates():
    sq16 = square_residues(16)
    table = {}
    for sign in (1, -1):
        for up, vp in ((1, 1), (1, 0), (0, 1)):
            residues = {
                (sign * 2 * u**4 + 25 * u*u*v*v + sign * 72 * v**4) % 16
                for u, v in product(range(16), repeat=2)
                if u % 2 == up and v % 2 == vp
            }
            assert not residues & sq16
            table[sign, up, vp] = residues
    assert table == {
        (1, 1, 1): {3, 11}, (-1, 1, 1): {7, 15},
        (1, 1, 0): {2, 6}, (-1, 1, 0): {2, 14},
        (1, 0, 1): {8, 12}, (-1, 0, 1): {8, 12},
    }
    sq49 = square_residues(49)
    admissible = [
        (u, v) for u, v in product(range(49), repeat=2)
        if (u % 7 or v % 7)
        and (7*u**4 - 50*u*u*v*v + 7*v**4) % 49 in sq49
    ]
    assert not admissible


def rhs(x):
    return x * (x + 9) * (x + 16)


def add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    xx, yy = q
    if x == xx and y == -yy:
        return None
    if p == q:
        slope = (3*x*x + 50*x + 144) / (2*y)
    else:
        slope = (yy - y) / (xx - x)
    out_x = slope*slope - 25 - x - xx
    return out_x, -y + slope*(x - out_x)


def check_points_and_reduction():
    points = {None} | {(F(x), F(y)) for x, y in [
        (0, 0), (-9, 0), (-16, 0), (12, 84), (12, -84),
        (-12, 12), (-12, -12),
    ]}
    assert len(points) == 8
    for p in points - {None}:
        x, y = p
        assert y*y == rhs(x)
    assert all(add(p, q) in points for p, q in product(points, repeat=2))
    assert add((F(12), F(84)), (F(12), F(84))) == (0, 0)
    assert add((F(-12), F(12)), (F(-12), F(12))) == (0, 0)
    discriminant = 16 * 144**2 * 49
    for prime in (5, 11):
        assert gcd(prime, discriminant) == 1
        count = 1 + sum(
            (y*y-rhs(x)) % prime == 0
            for x, y in product(range(prime), repeat=2)
        )
        assert count == 8


if __name__ == "__main__":
    check_local_certificates()
    check_points_and_reduction()
    print("Passed complete mod-16 and mod-49 descent certificates, both good-reduction counts, and the exact eight-point torsion table.")
