#!/usr/bin/env python3
"""Targeted actual-circle scan for shared-prime cancellation.

This is exploratory finite data, not a proof.  It enumerates lattice points
in circles with coordinates bounded by ``BOUND``, tests consecutive windows
against the intrinsic endpoint coefficient

    C_arc = N^(1/4) * angular_width,

and removes the actual common Gaussian gcd of each window.  For the resulting
primitive tuple, a split prime p with norm exponent e has allocations
``a_i=v_pi(z_i)``.  The raw pair source weight is e log p and the reduced
pair direction has |a_i-a_j| log p; their exact cancellation deficit is
``(e-|a_i-a_j|) log p``.  The script reports this deficit and the angular-gap
spread for the best windows.

No endpoint family is asserted.  Floating point is used only to rank and
filter candidate angles; all Gaussian gcds, divisions, and valuations are
exact integers.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import atan2, gcd, isqrt, log, pi


Gaussian = tuple[int, int]


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def norm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def ggcd(z: Gaussian, w: Gaussian) -> Gaussian:
    while w != (0, 0):
        qnum = mul(z, conj(w))
        den = norm(w)
        q = tuple((2 * x + den) // (2 * den) for x in qnum)
        z, w = w, (z[0] - mul(q, w)[0], z[1] - mul(q, w)[1])
    return z


def exact_div(z: Gaussian, w: Gaussian) -> Gaussian:
    qnum = mul(z, conj(w))
    den = norm(w)
    assert qnum[0] % den == 0 and qnum[1] % den == 0
    return qnum[0] // den, qnum[1] // den


def factor_integer(n: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def gaussian_prime(p: int) -> Gaussian:
    if p == 2:
        return (1, 1)
    assert p % 4 == 1
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b and b * b == b2:
            return a, b
    raise AssertionError(p)


def valuation(z: Gaussian, prime: Gaussian) -> int:
    p = norm(prime)
    value = 0
    while True:
        qnum = mul(z, conj(prime))
        if qnum[0] % p or qnum[1] % p:
            return value
        z = qnum[0] // p, qnum[1] // p
        value += 1


def primitive_tuple(points: tuple[Gaussian, ...]) -> tuple[tuple[Gaussian, ...], Gaussian]:
    common = points[0]
    for z in points[1:]:
        common = ggcd(common, z)
    return tuple(exact_div(z, common) for z in points), common


def circle_points(bound: int) -> dict[int, list[Gaussian]]:
    circles: dict[int, list[Gaussian]] = defaultdict(list)
    for x in range(-bound, bound + 1):
        for y in range(-bound, bound + 1):
            if x or y:
                circles[x * x + y * y].append((x, y))
    return circles


def profile(points: tuple[Gaussian, ...]) -> tuple[float, float, list[tuple[int, int, int]]]:
    """Return average cancellation fraction, gap ratio, and layer data."""
    n = norm(points[0])
    factors = factor_integer(n)
    split = [(p, e, gaussian_prime(p)) for p, e in factors.items()
             if p % 4 == 1]
    pair_raw = pair_reduced = 0.0
    for i, j in combinations(range(len(points)), 2):
        raw = reduced = 0.0
        for p, e, prime in split:
            ai, aj = valuation(points[i], prime), valuation(points[j], prime)
            assert 0 <= ai <= e and 0 <= aj <= e
            raw += e * log(p)
            reduced += abs(ai - aj) * log(p)
        pair_raw += raw
        pair_reduced += reduced
    cancellation = 0.0 if pair_raw == 0 else 1.0 - pair_reduced / pair_raw

    angles = [atan2(y, x) % (2 * pi) for x, y in points]
    gaps = [(angles[i + 1] - angles[i]) % (2 * pi)
            for i in range(len(angles) - 1)]
    gap_ratio = max(gaps) / min(gaps) if gaps else 1.0
    layers = []
    for p, e, prime in split:
        allocations = [valuation(z, prime) for z in points]
        for level in range(1, e + 1):
            layers.append((p, level, sum(a >= level for a in allocations)))
    return cancellation, gap_ratio, layers


def check_fixture() -> None:
    """Exact arithmetic regression for the displayed eight-point tuple."""
    points = ((70, -25), (71, -22), (73, -14), (74, -7),
              (74, 7), (73, 14), (71, 22), (70, 25))
    primitive, common = primitive_tuple(points)
    assert primitive == points and norm(common) == 1
    assert norm(points[0]) == 5525
    expected = {
        5: (2, [1, 0, 2, 2, 0, 0, 2, 1], 56, 30),
        13: (1, [1, 0, 1, 0, 1, 0, 1, 0], 28, 16),
        17: (1, [0, 1, 1, 0, 1, 0, 0, 1], 28, 16),
    }
    for p, (e, allocations, raw, reduced) in expected.items():
        actual = [valuation(z, gaussian_prime(p)) for z in points]
        assert actual == allocations
        assert raw == e * (len(points) * (len(points) - 1) // 2)
        assert reduced == sum(abs(a - b) for a, b in combinations(allocations, 2))
    tangents = []
    n = norm(points[0])
    for a, b in zip(points, points[1:]):
        cross = abs(a[0] * b[1] - a[1] * b[0])
        dot = a[0] * b[0] + a[1] * b[1]
        tangents.append(Fraction(cross, n + dot))
    assert tangents == [Fraction(1, 47), Fraction(1, 18), Fraction(1, 21),
                        Fraction(7, 74), Fraction(1, 21), Fraction(1, 18),
                        Fraction(1, 47)]
    assert max(tangents) / min(tangents) == Fraction(329, 74)


def scan(bound: int = 500, min_rows: int = 5, max_rows: int = 10,
         cutoff: float = 4.0) -> dict[int, tuple]:
    circles = circle_points(bound)
    best: dict[int, tuple] = {}
    candidates: dict[int, list[tuple]] = defaultdict(list)
    tested = 0
    for n, raw in circles.items():
        if len(raw) < min_rows:
            continue
        ordered = sorted(raw, key=lambda z: atan2(z[1], z[0]) % (2 * pi))
        doubled = ordered + ordered
        angles = [atan2(y, x) % (2 * pi) for x, y in ordered]
        doubled_angles = angles + [a + 2 * pi for a in angles]
        for rows in range(min_rows, min(max_rows, len(ordered)) + 1):
            for start in range(len(ordered)):
                window = tuple(doubled[start:start + rows])
                width = doubled_angles[start + rows - 1] - doubled_angles[start]
                if width <= 0:
                    continue
                raw_c = n ** 0.25 * width
                primitive, common = primitive_tuple(window)
                n_primitive = norm(primitive[0])
                intrinsic_c = n_primitive ** 0.25 * width
                if intrinsic_c > cutoff:
                    continue
                tested += 1
                cancellation, gap_ratio, layers = profile(primitive)
                key = (intrinsic_c, gap_ratio, -cancellation)
                old = best.get(rows)
                record = (key, n, common, primitive, raw_c, intrinsic_c,
                          cancellation, gap_ratio, layers)
                candidates[rows].append(record)
                if old is None or key < old[0]:
                    best[rows] = record
    print("scanned coordinate bound {}, {} candidate windows".format(bound, tested))
    for rows in sorted(best):
        key, n, common, points, raw_c, intrinsic_c, cancellation, gap_ratio, layers = best[rows]
        print("m={}: N={}, primitive_N={}, raw_C={:.6g}, intrinsic_C={:.6g}".format(
            rows, n, norm(points[0]), raw_c, intrinsic_c
        ))
        print("  common_gcd={}, cancellation_fraction={:.6f}, gap_ratio={:.6f}".format(
            common, cancellation, gap_ratio
        ))
        print("  points={}".format(points))
        print("  layer_counts={}".format(layers))
        fair = (2 ** (rows - 2) - 1) / (2 ** (rows - 1) - 1)
        most_cancelled = max(candidates[rows], key=lambda r: r[6])
        print("  fair_binary_cancellation={:.6f}; max_actual={:.6f} at C={:.6g}, gap_ratio={:.6f}".format(
            fair, most_cancelled[6], most_cancelled[5], most_cancelled[7]
        ))
        for threshold in (2.0, 4.0, 8.0, 12.0, 20.0):
            eligible = [r for r in candidates[rows] if r[5] <= threshold]
            if eligible:
                high = max(eligible, key=lambda r: r[6])
                low_gap = min(eligible, key=lambda r: r[7])
                print("  C<={:g}: max_cancel={:.6f}, min_gap_ratio={:.6f}, samples={}".format(
                    threshold, high[6], low_gap[7], len(eligible)
                ))
    return best


if __name__ == "__main__":
    check_fixture()
    print("PASS: exact eight-point fixture valuations, deficits, and gap tangents.")
    scan()
