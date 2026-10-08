#!/usr/bin/env python3
"""Exact checks for positive_pentagon_primitive_content.md."""

from fractions import Fraction
from itertools import combinations
from math import atan2, gcd, pi, prod, sqrt


def pentagon(a, b, c, d):
    assert 0 < a < b and 0 < c < d
    assert gcd(a, b) == gcd(c, d) == 1
    e = a * d + b * c - b * d
    assert e > 0
    H, X, Y = gcd(b, d), gcd(a, d - c), gcd(c, b - a)
    assert gcd(H, X) == gcd(H, Y) == gcd(X, Y) == 1
    assert gcd(b * (d - c), a * d) == H * X
    assert gcd(e, a * c) == X * Y
    assert gcd(d * (b - a), b * c) == H * Y
    assert e % (H * X * Y) == 0
    F = e // (H * X * Y)
    triples = [
        (a, b - a, b),
        (c, d - c, d),
        (b * (d - c) // (H * X), Y * F, a * d // (H * X)),
        (H * F, (b - a) * (d - c) // (X * Y), a * c // (X * Y)),
        (d * (b - a) // (H * Y), X * F, b * c // (H * Y)),
    ]
    u = []
    for n, r, den in triples:
        assert n > 0 and r > 0 and n + r == den
        assert gcd(n, den) == gcd(n, r) == gcd(r, den) == 1
        u.append(Fraction(n, den))
    for i in range(5):
        assert u[(i - 1) % 5] * u[(i + 1) % 5] == 1 - u[i]
    assert gcd((1 - u[2]).numerator, u[3].numerator,
               (1 - u[4]).numerator) == F
    return u, (H, X, Y), F


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def sub(z, w):
    return z[0] - w[0], z[1] - w[1]


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gaussian_gcd(z, w):
    while w != (0, 0):
        den = norm(w)
        # Exact nearest integers; no floating point in Euclidean division.
        real = z[0] * w[0] + z[1] * w[1]
        imag = z[1] * w[0] - z[0] * w[1]
        q = ((2 * real + den) // (2 * den),
             (2 * imag + den) // (2 * den))
        z, w = w, sub(z, mul(q, w))
    return z


def group_norm_gcd(points):
    return norm(group_gcd(points))


def group_gcd(points):
    g = (0, 0)
    for z in points:
        g = gaussian_gcd(g, z)
    return g


def gaussian_div(z, w):
    a, b = mul(z, (w[0], -w[1]))
    d = norm(w)
    assert a % d == b % d == 0
    return a // d, b // d


def valuation(z, prime):
    out = 0
    while z != (0, 0):
        a, b = mul(z, (prime[0], -prime[1]))
        d = norm(prime)
        if a % d or b % d:
            return out
        z = a // d, b // d
        out += 1


def integer_valuation(n, p):
    out = 0
    while n and n % p == 0:
        n //= p
        out += 1
    return out


def integer_content(vectors):
    out = 0
    for x, y in vectors:
        out = gcd(out, gcd(abs(x), abs(y)))
    return out


def product_ratio(z, pairs_n, pairs_d):
    A = mul(*(sub(z[i], z[j]) for i, j in pairs_n))
    B = mul(*(sub(z[i], z[j]) for i, j in pairs_d))
    assert A[0] * B[1] == A[1] * B[0]
    q = abs(Fraction(A[0], B[0]) if B[0] else Fraction(A[1], B[1]))
    # The parity-based primitive numerator bound, without square roots.
    assert 4 * q.numerator ** 2 <= norm(A)
    return q


def triangle(z, i, j, k):
    u, v = sub(z[j], z[i]), sub(z[k], z[i])
    return u[0] * v[1] - u[1] * v[0]


def circle_fixture(z, expected):
    N = norm(z[0])
    assert all(norm(w) == N for w in z)
    assert len(set(z)) == 5 and group_norm_gcd(z) == 1 and N % 2 == 1
    u = [product_ratio(z,
                      [(i, (i + 3) % 5), ((i + 1) % 5, (i + 2) % 5)],
                      [(i, (i + 2) % 5), ((i + 1) % 5, (i + 3) % 5)])
         for i in range(5)]
    found, HXY, F = pentagon(u[0].numerator, u[0].denominator,
                            u[1].numerator, u[1].denominator)
    assert found == u
    AO = group_norm_gcd([z[0], z[4]])
    AI = group_norm_gcd(z[1:4])
    NO, NI = N // AO, N // AI
    K = N // gcd(N, NO * NI)
    assert F % K == 0 and AI % K == 0
    Ds = {I: abs(triangle(z, *I)) for I in combinations(range(5), 3)}
    assert all(D > 0 and D % 2 == 0 for D in Ds.values())
    g = gcd(*Ds.values())
    g0 = gcd(*(D for I, D in Ds.items() if 4 not in I))
    g1 = gcd(*(D for I, D in Ds.items() if 0 not in I))
    assert (g * Ds[1, 2, 3]) % (g0 * g1) == 0
    Q = g * Ds[1, 2, 3] // (g0 * g1)
    assert 2 * K <= 2 * AI <= Ds[1, 2, 3]
    assert (N, u[0], u[1], HXY, F, K, Q) == expected
    for a, b, ratio in [(2, 3, 1 - u[2]), (1, 3, u[3]),
                        (1, 2, 1 - u[4])]:
        direct = product_ratio(z, [(0, 4), (a, b)], [(0, b), (a, 4)])
        assert ratio == direct and direct.numerator % K == 0
        assert 4 * F * F <= norm(sub(z[0], z[4])) * norm(sub(z[a], z[b]))
    # This angle calculation is diagnostic only; exact identities above do
    # not depend on floating point. These fixtures are outside C=1/2.
    angles = [atan2(w[1], w[0]) for w in z]
    steps = [(angles[i + 1] - angles[i] + pi) % (2 * pi) - pi
             for i in range(4)]
    assert all(s > 0 for s in steps) or all(s < 0 for s in steps)
    Delta = abs(sum(steps))
    assert Delta < pi
    C = Delta * sqrt(sqrt(N))
    assert C > 0.5
    # A C=1/2 endpoint arc would give |z_4-z_0| <= sqrt(R)/2,
    # hence 16 |z_4-z_0|^4 <= N.  Exclude it with exact integers.
    assert 16 * norm(sub(z[4], z[0])) ** 2 > N
    print(f"fixture N={N}: F={F}, K={K}, Q={Q}, diagnostic C={C:.8f}")


def source_gap_fixture(z, prime_powers, expected):
    """Check the exact gap valuation on a literal Gaussian circle tuple."""
    N = norm(z[0])
    assert all(norm(w) == N for w in z)
    assert len(set(z)) == 5 and group_norm_gcd(z) == 1 and N % 2 == 1
    assert N == prod(norm(p) ** e for p, e in prime_powers)
    orientations = [triangle(z, *indices) for indices in combinations(range(5), 3)]
    assert all(v > 0 for v in orientations) or all(v < 0 for v in orientations)

    # The given order is cyclic.  Reconstruct F from the actual chord ratios.
    ratios = [product_ratio(
        z,
        [(i, (i + 3) % 5), ((i + 1) % 5, (i + 2) % 5)],
        [(i, (i + 2) % 5), ((i + 1) % 5, (i + 3) % 5)])
        for i in range(5)]
    found, _, F = pentagon(ratios[0].numerator, ratios[0].denominator,
                           ratios[1].numerator, ratios[1].denominator)
    assert found == ratios

    G_O, G_I = group_gcd([z[0], z[4]]), group_gcd(z[1:4])
    x0, x4 = gaussian_div(z[0], G_O), gaussian_div(z[4], G_O)
    ys = [gaussian_div(z[i], G_I) for i in range(1, 4)]
    c_O = integer_content([sub(x4, x0)])
    c_I = integer_content([sub(ys[1], ys[0]), sub(ys[2], ys[0])])

    gaps = []
    units_by_row = []
    row_levels = [[] for _ in z]
    for prime, exponent in prime_powers:
        levels = [valuation(w, prime) for w in z]
        assert min(levels) == 0 and max(levels) == exponent
        for i, level in enumerate(levels):
            row_levels[i].append(level)
        outer = [levels[0], levels[4]]
        inner = levels[1:4]
        h = max(0, min(outer) - max(inner), min(inner) - max(outer))
        if not h:
            continue
        p = norm(prime)
        assert integer_valuation(F, p) == (
            h + integer_valuation(c_O, p) + integer_valuation(c_I, p))
        # Check the identity directly on all three reduced numerators too.
        for q in (1 - ratios[2], ratios[3], 1 - ratios[4]):
            assert integer_valuation(q.numerator, p) >= h
        gaps.append((p, h, tuple(outer), tuple(inner),
                     integer_valuation(c_O, p), integer_valuation(c_I, p),
                     integer_valuation(F, p)))

    # Recover the literal row units after removing every chosen source factor.
    for w, levels in zip(z, row_levels):
        source = (1, 0)
        for (prime, exponent), level in zip(prime_powers, levels):
            for _ in range(level):
                source = mul(source, prime)
            for _ in range(exponent - level):
                source = mul(source, (prime[0], -prime[1]))
        units_by_row.append(gaussian_div(w, source))
    assert all(u in ((1, 0), (0, 1), (-1, 0), (0, -1))
               for u in units_by_row)
    assert len(set(units_by_row)) >= 2
    assert (N, F, c_O, c_I, tuple(gaps)) == expected
    print(f"gap fixture N={N}: F={F}, cO={c_O}, cI={c_I}, gaps={gaps}")


def main():
    count = 0
    for b in range(2, 61):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            for d in range(2, 61):
                for c in range(1, d):
                    if gcd(c, d) == 1 and a * d + b * c > b * d:
                        pentagon(a, b, c, d)
                        count += 1
    assert count == 605550
    circle_fixture([(4, -33), (9, -32), (12, -31), (23, -24), (24, -23)],
                   (1105, Fraction(1, 2), Fraction(51, 52), (2, 1, 1), 25, 5, 5))
    circle_fixture([(-110, -35), (-109, -38), (-98, -61), (-94, -67), (-86, -77)],
                   (13325, Fraction(40, 41), Fraction(1, 2), (1, 1, 1), 39, 13, 13))
    # Literal Gaussian tuples, generated from the displayed split primes.
    # They include changing row units, both gap orientations, and nested powers.
    source_gap_fixture(
        [(4, -33), (9, -32), (12, -31), (23, -24), (24, -23)],
        [((2, 1), 1), ((3, 2), 1), ((4, 1), 1)],
        (1105, 25, 10, 1, ((5, 1, (1, 1), (0, 0, 0), 1, 0, 2),)))
    source_gap_fixture(
        [(-524, -7), (-488, -191), (-208, -481), (377, 364), (320, 415)],
        [((2, 1), 3), ((3, 2), 3)],
        (274625, 5275, 422, 5, ((5, 1, (3, 1), (0, 0, 0), 0, 1, 2),)))
    source_gap_fixture(
        [(-268, -1), (-257, -76), (-127, -236), (188, -191), (191, -188)],
        [((2, 1), 2), ((3, 2), 2), ((4, 1), 1)],
        (71825, 2125, 17, 5, ((5, 2, (0, 0), (2, 2, 2), 0, 1, 3),)))
    source_gap_fixture(
        [(-145, -10), (-143, -26), (-142, -31), (-130, -65), (-110, -95)],
        [((2, 1), 3), ((3, 2), 2)],
        (21125, 13, 1, 1, ((13, 1, (2, 2), (1, 0, 1), 0, 0, 1),)))
    for m in range(1, 1001):
        u, HXY, F = pentagon(3 * m + 1, 4 * m + 1, 3 * m + 1, 4 * m + 1)
        assert HXY == (4 * m + 1, 1, 1) and F == 2 * m + 1
        assert all(Fraction(1, 5) < x < Fraction(19, 20) for x in u)
    print(f"PASS: {count} reduced pentagons, exact Gaussian fixtures, "
          "gap-content valuations, primitive numerator bounds, and "
          "1000 shape-only examples.")


if __name__ == "__main__":
    main()
