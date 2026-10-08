"""Exact checks for the general near-real Gaussian-divisor reduction."""

from fractions import Fraction
from itertools import combinations, product
from math import gcd, prod, isqrt

from check_gaussian_reflection_replacement import (
    gmul, gconj, gnorm, gcd_many,
)


def power(z, exponent):
    out = (1, 0)
    for _ in range(exponent):
        out = gmul(out, z)
    return out


def tuple_from_columns(blocks, widths, rows):
    zs = []
    for row in rows:
        z = (1, 0)
        for block, width, exponent in zip(blocks, widths, row):
            z = gmul(z, power(block, exponent))
            z = gmul(z, power(gconj(block), width - exponent))
        zs.append(z)
    return zs


def pair_factor(blocks, a, b):
    out = (1, 0)
    for block, u, v in zip(blocks, a, b):
        out = gmul(out, power(block if u > v else gconj(block), abs(u-v)))
    return out


def associate_conjugacy_class(z):
    candidates = []
    for w in (z, gconj(z)):
        for _ in range(4):
            candidates.append(w)
            w = (-w[1], w[0])
    return min(candidates)


def valuation_two(value):
    value = abs(value)
    assert value != 0
    return (value & -value).bit_length()-1


def check_quotients(blocks, widths, rows):
    zs = tuple_from_columns(blocks, widths, rows)
    n = prod(gnorm(b) ** s for b, s in zip(blocks, widths))
    assert all(gnorm(z) == n for z in zs)
    assert gnorm(gcd_many(zs)) == 1
    factors = []
    for a, b in combinations(range(len(rows)), 2):
        factor = pair_factor(blocks, rows[a], rows[b])
        p = gnorm(factor)
        assert n % p == 0
        assert factor[1] != 0
        assert gnorm(gcd_many((factor, gconj(factor)))) == 1
        assert gmul(zs[a], gconj(factor)) == gmul(zs[b], factor)
        factors.append(factor)
    # The weighted Plotkin inequality, with logarithms removed.
    assert prod(gnorm(a) for a in factors) ** 4 <= n ** (len(rows) ** 2)
    for ia, ib, ic in combinations(range(len(rows)), 3):
        a = pair_factor(blocks, rows[ia], rows[ib])
        b = pair_factor(blocks, rows[ib], rows[ic])
        c = pair_factor(blocks, rows[ia], rows[ic])
        q = gnorm(gcd_many((a, gconj(b))))
        assert q % 2 == 1 and n % q == 0
        assert gmul(a, b) == (q*c[0], q*c[1])
        assert gnorm(a)*gnorm(b) == q*q*gnorm(c)
        ka, kb, kc = [valuation_two(z[1]) for z in (a, b, c)]
        assert kc == min(ka, kb) if ka != kb else kc > ka
    return zs, n, factors


def check_layer_cake():
    cases = 0
    for m in range(2, 8):
        for s in range(1, 4):
            for values in product(range(s+1), repeat=m):
                pair_sum = sum(abs(a-b) for a, b in combinations(values, 2))
                layers = [sum(e >= t for e in values) for t in range(1, s+1)]
                assert pair_sum == sum(k * (m-k) for k in layers)
                assert 4 * pair_sum <= s * m*m
                cases += 1
    return cases


def check_pell(index):
    u, v = 1, 0
    for _ in range(index):
        u, v = 9*u + 20*v, 4*u + 9*v
    assert u*u - 5*v*v == 1 and index % 10 == 1
    blocks = [(-1, 2), (u+v, 2*v), (2*v, u-v), (2*u+8*v, 3*u+v)]
    rows = [(1, 1, 1, 1), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    zs, n, factors = check_quotients(blocks, [1]*4, rows)

    # The common minus sign of the published family does not affect angles.
    ordered = [zs[i] for i in (1, 3, 0, 2)]
    relative = [gmul(z, gconj(ordered[0])) for z in ordered]
    assert relative[0] == (n, 0)
    assert all(x > 0 and y >= 0 for x, y in relative)
    assert all(a[0]*b[1] - a[1]*b[0] > 0
               for a, b in combinations(relative, 2))
    x, y = relative[-1]
    tangent_half = Fraction(y, n+x)
    # Delta=2 atan(tangent_half)<2*tangent_half proves C<2 exactly.
    assert tangent_half > 0 and tangent_half**4 * n < 1

    m = len(zs)
    assert all(gnorm(a)**2 >= n for a in factors)
    # The individual angular estimate gives norm injectivity for EVERY pair,
    # including those outside the near-critical Plotkin interval.
    assert len({gnorm(a) for a in factors}) == m*(m-1)//2
    assert all(abs(a[0]) == isqrt(gnorm(a)) for a in factors)
    assert all(a[1]**4 * n <= gnorm(a)**2 for a in factors)
    lifted = set()
    for (z1, z2), factor in zip(combinations(zs, 2), factors):
        p = gnorm(factor)
        g = n // p
        squared = gmul(factor, factor)
        x, y = (g*squared[0], g*squared[1])
        assert (x, y) == gmul(z1, gconj(z2))
        assert x*x+y*y == n*n and x > 0 and y != 0
        assert gcd(x, y) == g
        h = n-x
        assert h == 2*g*factor[1]**2
        assert h*(2*n-h) == y*y and gcd(h, 2*n-h) == 2*g
        assert h*h <= 4*n  # C=2: h <= 2 sqrt(N).
        lifted.add((x, abs(y)))
    assert len(lifted) == m*(m-1)//2
    good = [a for a in factors if gnorm(a)**(2*m) <= n**(m+2)]
    classes = {associate_conjugacy_class(a) for a in good}
    assert len(classes) == len(good)
    assert len({gnorm(a) for a in good}) == len(good)
    assert len(good) >= (m*(m-2)+3)//4
    # For C=2, the lemma's imaginary-part upper bound has this integral form.
    assert all(1 <= abs(a[1]) and abs(a[1])**(2*m) <= n for a in good)
    assert all(abs(a[0]) == isqrt(gnorm(a)) for a in good)
    return len(good)


def check_norm_uniqueness():
    for d in range(2, 2001):
        h = isqrt(isqrt(d))
        assert h**4 <= d
        representations = []
        for b in range(1, h+1):
            a = isqrt(d-b*b)
            if a*a + b*b == d:
                representations.append((a, b))
                assert a == isqrt(d)
        assert len(representations) <= 1


def check_binary_minimum():
    for m in range(2, 65):
        expected = 0
        cells = 2
        while cells < m:
            q, r = divmod(m, cells)
            expected += cells*q*(q-1)//2+r*q
            cells *= 2
        actual = sum(valuation_two(j-i)
                     for i, j in combinations(range(m), 2))
        assert actual == expected
        if m & (m-1) == 0:
            t = m.bit_length()-1
            assert 2*expected == m*m-m-m*t
        assert m*m//4 >= (m*(m-2)+3)//4


def main():
    cases = check_layer_cake()
    # This tuple has genuine three- and four-level allocations. Its angles
    # are unrestricted; here we test only the general quotient/Plotkin identities.
    check_quotients([(2, 1), (3, 2)], [2, 3],
                   [(0, 0), (2, 1), (1, 3), (0, 2)])
    counts = [check_pell(index) for index in (1, 11, 21)]
    check_norm_uniqueness()
    check_binary_minimum()
    print(f"PASS: {cases} exact layer-cake cases and a nested-prime Gaussian tuple.")
    print(f"PASS: three actual C<2 Pell clusters, with good divisor-class counts {counts}.")
    print("PASS: norm uniqueness through 2000; all six pair norms and individual angular bounds in each cluster.")
    print("PASS: exact autocorrelation lifts, endpoint equations, coordinate gcds, and distinct lifted points.")
    print("PASS: exact Gaussian triangle cancellation and binary residue valuations in all four tuples.")
    print("PASS: binary occupancy minima, power-of-two formula, and selected-edge quota through M=64.")


if __name__ == "__main__":
    main()
