#!/usr/bin/env python3
"""Finite scan of primitive quartet Gaussian contents on small complete circles.

For every circle x^2+y^2=N <= 2000, and every consecutive ordered cluster of
5--8 points, this checks each four-point subset.  Each full cluster is
normalized once by its actual Gaussian gcd; quartets retain those coordinates.
For each quartet we compute

    U=(z1-z0)(z3-z2), V=(z2-z0)(z3-z1), gamma=gcd_Z[i](U,V),

and verify Norm(gamma)=Pi/s^2 from the moving-quadrilateral identity.  The
reported gcd is across quartet gamma norms, while the minimum is an individual
quartet minimum; they need not agree.  This is finite evidence, not a theorem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
from math import atan2, gcd, isqrt

from check_four_vertical_lcm_search import (exact_div, gcd_gaussian, mul,
                                             norm)

Gaussian = tuple[int, int]
Point = Gaussian


def sub(a: Gaussian, b: Gaussian) -> Gaussian:
    return (a[0] - b[0], a[1] - b[1])


def determinant(a: Point, b: Point, c: Point) -> int:
    u, v = sub(b, a), sub(c, a)
    return u[0] * v[1] - u[1] * v[0]


def circle_points(n: int) -> list[Point]:
    out: set[Point] = set()
    for x in range(-isqrt(n), isqrt(n) + 1):
        y2 = n - x * x
        if y2 < 0:
            continue
        y = isqrt(y2)
        if y * y == y2:
            out.add((x, y))
            out.add((x, -y))
    return sorted(out, key=lambda z: atan2(z[1], z[0]))


def tuple_gcd(points: tuple[Point, ...]) -> Gaussian:
    g: Gaussian = points[0]
    for z in points[1:]:
        g = gcd_gaussian(g, z)
    return g


def primitive_tuple(points: tuple[Point, ...]) -> tuple[Point, ...]:
    g = tuple_gcd(points)
    return tuple(exact_div(z, g) for z in points)


def quartet_gamma(points: tuple[Point, Point, Point, Point]) -> tuple[int, Gaussian, int, int]:
    p0, p1, p2, p3 = points
    a, b, c, e = (determinant(p0, p1, p2), determinant(p0, p1, p3),
                  determinant(p0, p2, p3), determinant(p1, p2, p3))
    assert min(a, b, c, e) > 0 and e == a - b + c
    pi = a * b * c * e
    d01 = norm(sub(p1, p0))
    d02 = norm(sub(p2, p0))
    m = gcd(b * d02, c * d01)
    assert a * b * c % m == 0
    aa, bb, s = b * d02 // m, c * d01 // m, a * b * c // m
    assert aa > bb > 0 and gcd(aa, bb) == 1
    u = mul(sub(p1, p0), sub(p3, p2))
    v = mul(sub(p2, p0), sub(p3, p1))
    gamma = gcd_gaussian(u, v)
    assert gamma != (0, 0)
    assert pi % (s * s) == 0
    assert norm(gamma) == pi // (s * s)
    return norm(gamma), gamma, s, pi


def main() -> None:
    by_m: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    examples: dict[int, tuple[int, tuple[Point, ...], int, Gaussian]] = {}
    counterexamples: dict[int, tuple[int, int, tuple[Point, ...], int, Gaussian]] = {}
    primitive_quartets = 0
    clusters = 0

    for n in range(1, 2001):
        ordered = circle_points(n)
        if len(ordered) < 5:
            continue
        doubled = ordered + ordered
        for m in range(5, min(8, len(ordered)) + 1):
            for start in range(len(ordered)):
                cluster = tuple(doubled[start:start + m])
                clusters += 1
                # Normalize the complete source cluster once.  Per-quartet
                # normalization would erase the common Gaussian content.
                points = primitive_tuple(cluster)
                assert norm(tuple_gcd(points)) == 1
                values: list[int] = []
                gammas: list[Gaussian] = []
                inds_list = list(combinations(range(m), 4))
                for inds in inds_list:
                    quartet = tuple(points[i] for i in inds)
                    value, gamma, _, _ = quartet_gamma(quartet)
                    primitive_quartets += 1
                    values.append(value)
                    gammas.append(gamma)
                    old = examples.get(m)
                    if old is None or value < old[0]:
                        examples[m] = (value, quartet, n, gamma)
                g_int = gcd(*values)
                g_gauss: Gaussian = gammas[0]
                for gamma in gammas[1:]:
                    g_gauss = gcd_gaussian(g_gauss, gamma)
                # Circle-specific conjugation structure is expected to make
                # the Gaussian common-gcd norm equal the ordinary gcd of the
                # quartet norms.  Keep this as an explicit falsification test.
                assert norm(g_gauss) == g_int, (n, points, values, gammas,
                                                norm(g_gauss), g_int)
                by_m[m].append((min(values), g_int, norm(g_gauss)))
                # A cluster-level counterexample to identifying its minimum
                # with its common gcd.
                if min(values) != g_int and m not in counterexamples:
                    j = values.index(min(values))
                    counterexamples[m] = (min(values), g_int,
                                          tuple(points[i] for i in inds_list[j]),
                                          n, gammas[j])

    print(f"scanned {clusters:,} primitive source clusters and {primitive_quartets:,} quartet occurrences")
    for m in sorted(by_m):
        records = by_m[m]
        mins, gcds, ggcds = zip(*records)
        print(f"m={m}: clusters={len(records):,} "
              f"min_range=({min(mins)},{max(mins)}) "
              f"cluster_gcd_range=({min(gcds)},{max(gcds)}) "
              f"gaussian_gcd_norm_range=({min(ggcds)},{max(ggcds)})")
        print(f"     hist cluster_gcd={Counter(gcds).most_common(8)}; "
              f"gaussian_gcd_norm={Counter(ggcds).most_common(8)}")
        value, points, n, gamma = examples[m]
        print(f"     min example: N={n}, points={points}, gamma={gamma}, norm={value}")
        if m in counterexamples:
            value, cluster_gcd, points, n, gamma = counterexamples[m]
            print(f"     min-vs-cluster-gcd example: N={n}, points={points}, "
                  f"gamma={gamma}, min={value}, cluster_gcd={cluster_gcd}")


def mandatory_pi_family() -> None:
    """Check the requested primitive family, retaining its source content."""
    pi: Gaussian = (2, 1)
    print("mandatory pi=2+i family:")
    for e in range(1, 5):
        z = (1, 0)
        for _ in range(e):
            z = mul(z, pi)
        raw = (z, (-z[1], z[0]), (-z[0], -z[1]), (z[1], -z[0]), (z[0], -z[1]))
        raw = tuple(sorted(raw, key=lambda q: atan2(q[1], q[0])))
        points = primitive_tuple(raw)
        # The gcd is a unit; primitive_tuple may rotate by that unit, which
        # preserves the source content and all gamma norms.
        assert norm(tuple_gcd(points)) == 1
        vals = []
        gammas = []
        for inds in combinations(range(5), 4):
            value, gamma, _, _ = quartet_gamma(tuple(points[i] for i in inds))
            vals.append(value)
            gammas.append(gamma)
        gg = gammas[0]
        for gamma in gammas[1:]:
            gg = gcd_gaussian(gg, gamma)
        assert norm(gg) == gcd(*vals)
        print(f"  e={e}, N={norm(z)}, min={min(vals)}, gcd_norm={gcd(*vals)}, "
              f"gaussian_gcd_norm={norm(gg)}, norms={sorted(set(vals))}")


if __name__ == "__main__":
    # Run the family first so a regression in its source-content behavior is
    # immediately visible.
    mandatory_pi_family()
    main()
