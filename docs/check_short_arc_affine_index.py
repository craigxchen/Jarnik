#!/usr/bin/env python3
"""Bounded affine-index/triangle checks on complete integer circles."""

from itertools import combinations
from math import comb, gcd, isqrt


def points(n):
    r = isqrt(n)
    return [(x, y) for x in range(-r, r + 1) for y in range(-r, r + 1)
            if x * x + y * y == n]


def det(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def affine_data(P):
    ds = [abs(det(*T)) for T in combinations(P, 3)]
    g = 0
    for d in ds:
        g = gcd(g, d)
    primitive = [d // g for d in ds]
    return g, max(primitive), primitive


def arc_windows(P, n, m, C):
    # Sort by angle, then inspect cyclic windows.  The chord bound is a
    # sufficient short-arc filter and avoids floating point angle tests.
    P = sorted(P, key=lambda z: __import__('math').atan2(z[1], z[0]))
    out = []
    for i in range(len(P)):
        W = [P[(i + j) % len(P)] for j in range(m)]
        if len(set(W)) < m:
            continue
        max_chord2 = max((x - y) ** 2 + (u - v) ** 2
                         for (x, u), (y, v) in combinations(W, 2))
        # C*sqrt(R) has squared length C^2*sqrt(N).
        if max_chord2 <= C * C * n ** 0.5:
            out.append(W)
    return out


if __name__ == "__main__":
    fixtures = [5, 25, 65, 85, 125, 325, 425]
    seen = 0
    for n in fixtures:
        P = points(n)
        print(f"N={n}: r2={len(P)}")
        for m in (5, 6):
            if len(P) < m:
                continue
            windows = arc_windows(P, n, m, 20.0)
            best = None
            for W in windows:
                g, H, primitive = affine_data(W)
                b = comb(m // 2, 3) + comb((m + 1) // 2, 3)
                product = 1
                for q in primitive:
                    product *= q
                assert product % (n ** b) == 0
                candidate = (H, g, W)
                best = candidate if best is None or candidate[:2] < best[:2] else best
                seen += 1
            if best:
                print(f"  m={m}: chord-filtered windows={len(windows)}, min(H,g)=({best[0]},{best[1]})")
    assert seen > 0
    print("PASS: triangle-content divisibility and bounded short-window scan")
