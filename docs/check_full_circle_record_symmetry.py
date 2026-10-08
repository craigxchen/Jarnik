#!/usr/bin/env python3
"""Exact finite checks for the full-circle record symmetry argument.

This is a bounded sanity check, not a proof of the unbounded record claim.
All phase and projective computations use ``Fraction``; there are no external
dependencies.
"""

from fractions import Fraction
from functools import cmp_to_key
from math import isqrt


def circle_points(n):
    bound = isqrt(n)
    return [(x, y) for x in range(-bound, bound + 1)
            for y in range(-bound, bound + 1) if x * x + y * y == n]


def phase(n, z0, z):
    """Return the normalized phase z/z0, as an integer pair over n."""
    x0, y0 = z0
    x, y = z
    return (x * x0 + y * y0, y * x0 - x * y0)


def cyclic_order(phases):
    def half(p):
        x, y = p
        return y > 0 or (y == 0 and x >= 0)

    def compare(a, b):
        ha, hb = half(a), half(b)
        if ha != hb:
            return -1 if ha else 1
        cross = a[0] * b[1] - a[1] * b[0]
        return (cross > 0) - (cross < 0)

    return sorted(phases, key=cmp_to_key(compare))


def half_angle(p, n):
    """Homogeneous half-angle coordinates h=(u,v), t=v/u."""
    x, y = p
    # For q=(x+iy)/n, q=(1+it)/(1-it), hence t=y/(n+x).
    if n + x == 0 and y == 0:
        return (Fraction(0), Fraction(1))
    return (Fraction(n + x), Fraction(y))


def mobius_from_triples(source, target):
    rows = []
    for (u, v), (U, V) in zip(source, target):
        rows.append([-V * u, -V * v, U * u, U * v])
    # RREF of the three homogeneous linear constraints, then take a free
    # variable as 1.  A projective map is unique when the labels are distinct.
    a = [list(row) + [Fraction(0)] for row in rows]
    pivots = []
    r = 0
    for c in range(4):
        pivot = next((i for i in range(r, 3) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [x / scale for x in a[r]]
        for i in range(3):
            if i != r and a[i][c]:
                scale = a[i][c]
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
    free = next(c for c in range(4) if c not in pivots)
    vec = [Fraction(0)] * 4
    vec[free] = Fraction(1)
    for i, c in enumerate(pivots):
        vec[c] = -a[i][free]
    return tuple(vec)


def apply(m, h):
    a, b, c, d = m
    u, v = h
    return (a * u + b * v, c * u + d * v)


def same_projective(a, b):
    return a[0] * b[1] == a[1] * b[0]


def determinant(m):
    a, b, c, d = m
    return a * d - b * c


def q_coordinates(h):
    u, v = h
    return (u * u - v * v, 2 * u * v, u * u + v * v)


def classify(source_h, target_h, m):
    """Check whether the induced phase map is lambda*q or lambda*conj(q)."""
    def normalized(h):
        x, y, d = q_coordinates(h)
        return (x / d, y / d)
    src = [normalized(h) for h in source_h]
    dst = [normalized(apply(m, h)) for h in source_h]
    lx, ly = src[0]
    wx, wy = dst[0]
    for conjugate in (False, True):
        # lambda = target/source, or target/conj(source), respectively.
        lam = ((wx * lx + wy * ly, wy * lx - wx * ly) if not conjugate
               else (wx * lx - wy * ly, wy * lx + wx * ly))
        def image(s):
            sx, sy = s
            if conjugate:
                sy = -sy
            return (lam[0] * sx - lam[1] * sy,
                    lam[0] * sy + lam[1] * sx)
        if all(image(s) == d for s, d in zip(src, dst)):
            return "anticonformal" if conjugate else "conformal"
    return "NONCONFORMAL"


def check_circle(n):
    points = circle_points(n)
    z0 = points[0]
    ordered = cyclic_order([phase(n, z0, z) for z in points])
    hs = [half_angle(p, n) for p in ordered]
    count = 0
    kinds = set()
    k = len(hs)
    for shift in range(k):
        for reverse in (False, True):
            target = [hs[(shift + (-i if reverse else i)) % k] for i in range(k)]
            m = mobius_from_triples(hs[:3], target[:3])
            if determinant(m) == 0:
                continue
            if all(same_projective(apply(m, s), t) for s, t in zip(hs, target)):
                a, b, c, d = m
                assert a*a+c*c == b*b+d*d and a*b+c*d == 0
                kinds.add(classify(hs, target, m))
                count += 1
    assert count == 8 and kinds == {"conformal", "anticonformal"}
    return len(points), count, kinds


def records_upto(limit):
    counts = []
    best = 0
    for n in range(1, limit + 1):
        r = len(circle_points(n))
        if r > best:
            counts.append((n, r))
            best = r
    return counts


if __name__ == "__main__":
    fixtures = [1, 5, 13, 25, 65, 85, 325]
    for n in fixtures:
        points, maps, kinds = check_circle(n)
        print(f"N={n}: r2={points}, surviving maps={maps}, types={sorted(kinds)}")
    rec = records_upto(325)
    print("strict records up to 325:", rec)
    for n, r in rec:
        assert r > max((len(circle_points(j)) for j in range(1, n)), default=0)
    print("These finite checks supplement, rather than prove, the unbounded record theorem.")
    print("PASS: exact full-circle symmetry and bounded record checks")
