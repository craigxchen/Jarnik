#!/usr/bin/env python3
"""Exact finite comparison of integer and Gaussian cluster contents.

Each ordered m=5..8 cluster is normalized once by its common Gaussian gcd.
For the resulting primitive cluster, let g be the gcd of all triangle
determinants and G the Gaussian gcd of all quartet gamma values.  The scan
tests candidate relations such as N(G)|g^2 and records complete-cluster
witnesses.  It is finite experimentation only.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
from math import atan2, gcd, prod

from check_moving_quadrilateral_gamma_content import (circle_points,
                                                       determinant,
                                                       primitive_tuple,
                                                       quartet_gamma,
                                                       tuple_gcd)
from check_four_vertical_lcm_search import gcd_gaussian, mul, norm

Point = tuple[int, int]
Gaussian = Point


def cluster_data(points: tuple[Point, ...]) -> tuple[int, Gaussian, int, list[int]]:
    N = norm(points[0])
    g = 0
    for i, j, k in combinations(range(len(points)), 3):
        g = gcd(g, abs(determinant(points[i], points[j], points[k])))
    gammas: list[Gaussian] = []
    norms: list[int] = []
    for inds in combinations(range(len(points)), 4):
        value, gamma, _, _ = quartet_gamma(tuple(points[i] for i in inds))
        gammas.append(gamma)
        norms.append(value)
    G = gammas[0]
    for gamma in gammas[1:]:
        G = gcd_gaussian(G, gamma)
    NG = norm(G)
    assert norm(tuple_gcd(points)) == 1
    # Global difference content and leave-one-out source content.
    D = tuple_gcd(tuple((points[i][0] - points[0][0],
                         points[i][1] - points[0][1])
                        for i in range(1, len(points))))
    ND = norm(D)
    H = [tuple_gcd(tuple(points[j] for j in range(len(points)) if j != i))
         for i in range(len(points))]
    E = prod(norm(h) for h in H)
    # These are the proposed global divisibilities, checked per full cluster.
    assert N % E == 0
    assert g % ND == 0
    assert (4 * E * g * g) % (NG * ND) == 0
    assert (2 * E * g * g) % NG == 0
    return g, G, NG, norms


def report(name: str, records: list[tuple[int, int, int, int, tuple[Point, ...]]]) -> None:
    # (circle norm, g, Norm(G), gcd(norm gamma), cluster)
    print(f"{name}: clusters={len(records):,}")
    print(f"  g range=({min(r[1] for r in records)},{max(r[1] for r in records)}), "
          f"NormG range=({min(r[2] for r in records)},{max(r[2] for r in records)})")
    print(f"  ratios g^2/NormG: "
          f"min={min(r[1]*r[1]/r[2] for r in records):g}, "
          f"max={max(r[1]*r[1]/r[2] for r in records):g}")
    for label, pred in (
        ("NormG|g^2", lambda r: r[1] * r[1] % r[2] == 0),
        ("NormG<=g^2", lambda r: r[2] <= r[1] * r[1]),
        ("gcd(norms)=NormG", lambda r: r[3] == r[2]),
    ):
        bad = next((r for r in records if not pred(r)), None)
        print(f"  {label}: {'PASS' if bad is None else 'FAIL'}")
        if bad is not None:
            print(f"    witness N={bad[0]}, g={bad[1]}, NormG={bad[2]}, "
                  f"gcd_norms={bad[3]}, cluster={bad[4]}")
    print(f"  common gcd histogram: {Counter(r[2] for r in records).most_common(10)}")


def pi_family() -> list[tuple[int, int, int, int, tuple[Point, ...]]]:
    pi: Gaussian = (2, 1)
    out = []
    z: Gaussian = (1, 0)
    for e in range(1, 5):
        z = mul(z, pi)
        raw = (z, (-z[1], z[0]), (-z[0], -z[1]), (z[1], -z[0]),
               (z[0], -z[1]))
        points = tuple(sorted(raw, key=lambda q: atan2(q[1], q[0])))
        points = primitive_tuple(points)
        g, G, NG, norms = cluster_data(points)
        D = tuple_gcd(tuple((points[i][0] - points[0][0],
                             points[i][1] - points[0][1])
                            for i in range(1, len(points))))
        H = [tuple_gcd(tuple(points[j] for j in range(5) if j != i))
             for i in range(5)]
        E = prod(norm(h) for h in H)
        out.append((norm(points[0]), g, NG, gcd(*norms), points))
        print(f"pi family e={e}: g={g}, NormD={norm(D)}, E={E}, "
              f"NormG={NG}, gcd_norms={gcd(*norms)}, cluster={points}")
    return out


def main() -> None:
    by_m: dict[int, list[tuple[int, int, int, int, tuple[Point, ...]]]] = defaultdict(list)
    doubled_clusters = 0
    for n in range(1, 2001):
        ordered = circle_points(n)
        if len(ordered) < 5:
            continue
        doubled = ordered + ordered
        for m in range(5, min(8, len(ordered)) + 1):
            for start in range(len(ordered)):
                doubled_clusters += 1
                points = primitive_tuple(tuple(doubled[start:start + m]))
                g, G, NG, norms = cluster_data(points)
                by_m[m].append((norm(points[0]), g, NG, gcd(*norms), points))
    print(f"complete-circle scan: {doubled_clusters:,} clusters")
    for m in sorted(by_m):
        report(f"m={m}", by_m[m])
    report("mandatory pi family", pi_family())


if __name__ == "__main__":
    main()
