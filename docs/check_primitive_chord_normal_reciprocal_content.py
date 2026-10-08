#!/usr/bin/env python3
"""Exact checks for the primitive-chord normal and reciprocal-content audit."""

from __future__ import annotations

from math import atan2, gcd, isqrt, sqrt


def det(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[1] - a[1] * b[0]


def dot(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[0] + a[1] * b[1]


def norm(a: tuple[int, int]) -> int:
    return dot(a, a)


def gaussian_mul(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def circle_points(n: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for x in range(-isqrt(n), isqrt(n) + 1):
        y2 = n - x * x
        y = isqrt(y2)
        if y * y != y2:
            continue
        out.append((x, y))
        if y:
            out.append((x, -y))
    return out


def check_frame(
    p0: tuple[int, int], p1: tuple[int, int], rows: list[tuple[int, int]]
) -> tuple[int, int, int, int, int]:
    n = norm(p0)
    assert n > 0 and norm(p1) == n
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    g = gcd(abs(dx), abs(dy))
    assert g > 0
    w = (dx // g, dy // g)
    assert gcd(abs(w[0]), abs(w[1])) == 1
    s = norm(w)
    gg = g * s
    assert gg % 2 == 0
    t0 = gg // 2
    b = det(w, p0)
    assert (2 * b) % s == 0
    t = 2 * b // s
    assert 4 * n == s * (g * g + t * t)
    assert (2 * p0[0], 2 * p0[1]) == gaussian_mul(w, (-g, t))
    assert (2 * p1[0], 2 * p1[1]) == gaussian_mul(w, (g, t))
    assert dot(w, p0) == -t0 and det(w, p0) == b

    if s > 1:
        assert gcd(abs(w[0]), s) == gcd(abs(w[1]), s) == 1
        k = w[1] * pow(w[0], -1, s) % s
        assert (k * k + 1) % s == 0
    else:
        k = 0

    normal = (-w[1], w[0])
    multiplier = (-w[1], -w[0])  # -i conjugate(w)
    d_direct = 0
    d_increments = gcd(abs(b), t0)
    checked = 0
    for z in rows:
        assert norm(z) == n
        h = (z[0] - p0[0], z[1] - p0[1])
        pp, qq = dot(w, h), det(w, h)
        assert (w[0] * pp - w[1] * qq) % s == 0
        assert (w[1] * pp + w[0] * qq) % s == 0
        if s > 1:
            assert (pp - k * qq) % s == 0
        assert (pp * pp + qq * qq - gg * pp + 2 * b * qq) == 0
        assert (2 * b * qq) % s == 0

        x, y = dot(normal, z), det(normal, z)
        assert (x, y) == (b + qq, t0 - pp)
        assert (x, y) == gaussian_mul(multiplier, z)
        assert x * x + y * y == s * n
        assert y * y == t0 * t0 - 2 * b * qq - qq * qq
        assert (y - t0) * (y + t0) == -qq * (2 * b + qq)
        assert pp * (gg - pp) == qq * (2 * b + qq)

        d_direct = gcd(d_direct, abs(x), abs(y))
        d_increments = gcd(d_increments, abs(pp), abs(qq))
        checked += 1
    assert p0 in rows and p1 in rows
    assert d_direct == d_increments
    assert gcd(abs(b), t0) * 2 == s * gcd(g, abs(t))
    return checked, s, b, t0, d_direct


def exhaustive_endpoint_frames() -> tuple[int, int]:
    cases = 0
    rows_checked = 0
    seen: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for u in range(-3, 4):
        for v in range(-3, 4):
            if (u, v) == (0, 0) or gcd(abs(u), abs(v)) != 1:
                continue
            for g in range(1, 4):
                for t in range(-5, 6):
                    left = (-g * u - t * v, -g * v + t * u)
                    if left[0] % 2 or left[1] % 2:
                        continue
                    p0 = (left[0] // 2, left[1] // 2)
                    p1 = (p0[0] + g * u, p0[1] + g * v)
                    if norm(p0) == 0 or (p0, p1) in seen:
                        continue
                    seen.add((p0, p1))
                    count, *_ = check_frame(p0, p1, circle_points(norm(p0)))
                    cases += 1
                    rows_checked += count
    return cases, rows_checked


def literal_short_arc() -> None:
    p0, interior, p1 = (7, 6), (6, 7), (2, 9)
    assert norm(p0) == norm(interior) == norm(p1) == 85
    arc = sqrt(85) * (atan2(p1[1], p1[0]) - atan2(p0[1], p0[0]))
    assert 0 < arc < 2 * 85**0.25
    count, s, b, tt, d = check_frame(p0, p1, [p0, interior, p1])
    assert count == 3 and (s, b, tt, d) == (34, -51, 17, 1)
    w = (-5, 3)
    h = (interior[0] - p0[0], interior[1] - p0[1])
    pp, qq = dot(w, h), det(w, h)
    assert (pp, qq) == (8, -2)
    c = qq * (2 * b + qq)
    assert c == 208 and tt * tt == 289 and gcd(c, tt * tt) == 1
    assert gcd(abs(b), tt) == 17


def main() -> None:
    cases, rows = exhaustive_endpoint_frames()
    literal_short_arc()
    print(
        f"PASS: {cases} equal-norm endpoint frames, {rows} literal circle rows, "
        "and one primitive three-row C<2 content-loss fixture"
    )


if __name__ == "__main__":
    main()
