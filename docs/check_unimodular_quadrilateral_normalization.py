#!/usr/bin/env python3
"""Finite checker for bounded-determinant quadrilateral normalization.

This checks the elementary reduction used in the bounded affine type claim:
translate P0 to zero, send p=P1-P0 to (g,0) with a GL2(Z) matrix, and
center the first coordinate of q=P2-P0 by a shear.  The remaining coordinates
are then bounded in terms of the common triple-determinant bound B.

It is a randomized exact checker, not a proof.  All arithmetic is integer.
"""

from __future__ import annotations

import math
import random
from itertools import combinations

Vec = tuple[int, int]
Mat = tuple[tuple[int, int], tuple[int, int]]


def det(u: Vec, v: Vec) -> int:
    return u[0] * v[1] - u[1] * v[0]


def mat_vec(m: Mat, v: Vec) -> Vec:
    return (m[0][0] * v[0] + m[0][1] * v[1],
            m[1][0] * v[0] + m[1][1] * v[1])


def mat_mul(a: Mat, b: Mat) -> Mat:
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))  # type: ignore[return-value]


def egcd(a: int, b: int) -> tuple[int, int, int]:
    """Return x,y,g with x*a+y*b=g and g=gcd(a,b)>=0."""
    if b == 0:
        return (1 if a >= 0 else -1, 0, abs(a))
    x, y, g = egcd(b, a % b)
    return y, x - (a // b) * y, g


def map_primitive_axis(p: Vec) -> Mat:
    """Return U in GL2(Z) with U*p=(g,0), g=gcd(abs(p))."""
    x, y = p
    a, b, g = egcd(x, y)
    # First row is (a,b), hence evaluates to g.  The second row completes it
    # to determinant one: (-y/g, x/g).
    u: Mat = ((a, b), (-y // g, x // g))
    assert abs(det(u[0], u[1])) == 1
    assert mat_vec(u, p) == (g, 0)
    return u


def inv_unimodular(u: Mat) -> Mat:
    d = det(u[0], u[1])
    assert abs(d) == 1
    return ((u[1][1] // d, -u[0][1] // d),
            (-u[1][0] // d, u[0][0] // d))


def centered_residue(b: int, s: int) -> int:
    """Residue b modulo positive s in [-s//2, s//2]."""
    r = b % s
    if 2 * r > s:
        r -= s
    return r


def normalize(points: tuple[Vec, Vec, Vec, Vec]) -> tuple[tuple[Vec, Vec, Vec, Vec], Mat]:
    p0, p1, p2, p3 = points
    shifted = tuple((x - p0[0], y - p0[1]) for x, y in points)
    p = shifted[1]
    q = shifted[2]
    u0 = map_primitive_axis(p)
    q1 = mat_vec(u0, q)
    g = math.gcd(abs(p[0]), abs(p[1]))
    s = q1[1]
    assert s != 0  # no three collinear implies det(p,q) != 0
    if s < 0:
        # Flip the second row; p remains (g,0), while q's second coordinate
        # becomes positive.  Determinant -1 is allowed in GL2(Z).
        u0 = (u0[0], (-u0[1][0], -u0[1][1]))
        q1 = mat_vec(u0, q)
    assert q1[1] > 0 and q1[1] == abs(det(p, q)) // g
    k = (q1[0] - centered_residue(q1[0], q1[1])) // q1[1]
    shear: Mat = ((1, -k), (0, 1))
    u = mat_mul(shear, u0)
    out = tuple(mat_vec(u, v) for v in shifted)
    assert out[0] == (0, 0) and out[1] == (g, 0)
    assert out[2][0] == centered_residue(out[2][0], out[2][1])
    return out, u


def triple_bound(points: tuple[Vec, Vec, Vec, Vec]) -> int:
    return max(abs(det((points[j][0] - points[i][0], points[j][1] - points[i][1]),
                       (points[k][0] - points[i][0], points[k][1] - points[i][1])))
               for i, j, k in combinations(range(4), 3))


def no_three_collinear(points: tuple[Vec, Vec, Vec, Vec]) -> bool:
    return all(det((points[j][0] - points[i][0], points[j][1] - points[i][1]),
                   (points[k][0] - points[i][0], points[k][1] - points[i][1])) != 0
               for i, j, k in combinations(range(4), 3))


def random_gl2(rng: random.Random, steps: int = 12) -> Mat:
    u: Mat = ((1, 0), (0, 1))
    moves = (((1, 1), (0, 1)), ((1, -1), (0, 1)),
             ((1, 0), (1, 1)), ((1, 0), (-1, 1)),
             ((0, 1), (1, 0)), ((0, -1), (1, 0)))
    for _ in range(steps):
        u = mat_mul(rng.choice(moves), u)
    return u


def main() -> None:
    rng = random.Random(0x5EED)
    checked = 0
    for bnd in range(1, 13):
        # Enumerating a small box supplies many shapes; random affine images
        # make the input coordinates large without changing determinants.
        candidates = [(x, y) for x in range(-bnd - 1, bnd + 2)
                      for y in range(-bnd - 1, bnd + 2)]
        for _ in range(5000):
            base = tuple(rng.sample(candidates, 4))
            if not no_three_collinear(base) or triple_bound(base) > bnd:
                continue
            t = random_gl2(rng, steps=30)
            z = (rng.randrange(-10**6, 10**6), rng.randrange(-10**6, 10**6))
            image = tuple((z[0] + v[0], z[1] + v[1]) for v in
                          (mat_vec(t, p) for p in base))
            normalized, u = normalize(image)
            ui = inv_unimodular(u)
            assert mat_mul(u, ui) == ((1, 0), (0, 1))
            assert mat_mul(ui, u) == ((1, 0), (0, 1))
            assert all(mat_vec(ui, normalized[i]) ==
                       (image[i][0] - image[0][0], image[i][1] - image[0][1])
                       for i in range(4))
            assert triple_bound(normalized) == triple_bound(base) <= bnd
            g = math.gcd(abs(image[1][0] - image[0][0]),
                         abs(image[1][1] - image[0][1]))
            assert g <= bnd  # g divides det(p,q), which is nonzero and <= B
            # A deliberately loose explicit bound following the two determinant
            # equations: |y_i| <= B/g and |x_i| <= B^2+B for i=2,3.
            assert all(abs(v[1]) <= bnd for v in normalized[2:])
            assert all(abs(v[0]) <= bnd * bnd + bnd for v in normalized[2:])
            checked += 1
    assert checked >= 100
    print(f"checked {checked} random translated/unimodular quadrilaterals")


if __name__ == "__main__":
    main()
