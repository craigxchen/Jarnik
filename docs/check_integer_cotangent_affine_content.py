"""Exact checks for integer-cotangent affine and Gaussian contents.

The proof is in integer_cotangent_affine_content.md.  These finite tests
exercise signs, parameter contents, the prime two, and the Pell family.
"""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm, prod
from random import Random

from check_integer_cotangent_normalization import edge_norm, edge_quotient
from check_least_radius_formula import (
    conj,
    exact_div,
    gcd_all,
    gcd_gaussian,
    mul,
    norm,
)


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def determinant(a, b):
    return a[0] * b[1] - a[1] * b[0]


def valuation(n, p):
    n = abs(n)
    answer = 0
    while n and n % p == 0:
        n //= p
        answer += 1
    return answer


def prime_divisors(n):
    n = abs(n)
    answer = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            answer.append(p)
            while n % p == 0:
                n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        answer.append(n)
    return answer


def universal_rows(xs, scale):
    hs = [(x, scale) for x in xs]
    bars = [conj(h) for h in hs]
    row0 = reduce(mul, bars, (1, 0))
    rows = [row0]
    for i, h in enumerate(hs):
        rows.append(mul(h, reduce(mul, (bars[j] for j in range(len(xs)) if j != i), (1, 0))))
    return rows, bars


def least_radius_squared(xs, scale):
    edges = list(xs)
    edges.extend(edge_quotient(x, y, scale) for x, y in combinations(xs, 2))
    return lcm(*(edge_norm(q, scale) for q in edges))


def check_data(xs, scale):
    assert scale > 0 and len(xs) >= 2 and len(set(xs)) == len(xs)
    for x, y in combinations(xs, 2):
        assert (x * y + scale * scale) % (x - y) == 0

    rows, bars = universal_rows(xs, scale)
    gaussian_content = gcd_all(rows)
    content_norm = norm(gaussian_content)
    primitive = [exact_div(row, gaussian_content) for row in rows]
    radius_squared = least_radius_squared(xs, scale)
    factor_norms = [x * x + scale * scale for x in xs]

    assert all(norm(z) == radius_squared for z in primitive)
    assert norm(gcd_all(primitive)) == 1
    assert content_norm * radius_squared == prod(factor_norms)

    triangle_determinants = []
    for a, b, c in combinations(range(len(primitive)), 3):
        triangle_determinants.append(
            abs(determinant(sub(primitive[b], primitive[a]), sub(primitive[c], primitive[a])))
        )
    triangle_content = reduce(gcd, triangle_determinants)
    assert triangle_content > 0
    assert (4 * scale**3) % triangle_content == 0
    assert scale % gcd(triangle_content, radius_squared) == 0

    raw_anchored = []
    for i, j in combinations(range(len(xs)), 2):
        actual = abs(determinant(sub(rows[i + 1], rows[0]), sub(rows[j + 1], rows[0])))
        predicted = (
            4
            * scale**3
            * abs(xs[i] - xs[j])
            * prod(factor_norms[q] for q in range(len(xs)) if q not in (i, j))
        )
        assert actual == predicted
        raw_anchored.append(actual)
    assert reduce(gcd, raw_anchored) == content_norm * triangle_content

    # The selected-pair Gaussian divisor (12), including its exact norm.
    for i, j in combinations(range(len(xs)), 2):
        pair_gcd = gcd_gaussian(bars[i], bars[j])
        delta = xs[i] - xs[j]
        d = gcd((xs[i] * xs[i] + scale * scale) // delta, xs[i], scale, delta)
        assert norm(pair_gcd) == abs(delta) * abs(d)
        divisor = mul((0, 2 * scale), pair_gcd)
        divisor = mul(
            divisor,
            reduce(mul, (bars[q] for q in range(len(xs)) if q not in (i, j)), (1, 0)),
        )
        exact_div(divisor, gaussian_content)
        assert norm(divisor) % content_norm == 0

    # Check the sharper exponent distinction prime by prime.
    for p in prime_divisors(triangle_content):
        exponent = valuation(triangle_content, p)
        if p == 2:
            assert exponent <= 3 * valuation(scale, 2) + 2
        elif radius_squared % p == 0:
            assert exponent <= valuation(scale, p)
        else:
            assert exponent <= 3 * valuation(scale, p)

    return triangle_content, content_norm, radius_squared


def cleared_clique(base_xs, base_scale=1):
    multiplier = 1
    for x, y in combinations(base_xs, 2):
        multiplier = lcm(
            multiplier,
            Fraction(x * y + base_scale * base_scale, x - y).denominator,
        )
    return tuple(multiplier * x for x in base_xs), multiplier * base_scale


def main():
    checks = 0
    maximum_triangle_content = 0

    # Direct small integral cliques.
    for scale in range(1, 8):
        values = range(scale, 25)
        compatible = {
            (x, y): (x * y + scale * scale) % (x - y) == 0
            for x, y in combinations(values, 2)
        }
        for size in (2, 3, 4):
            for xs in combinations(values, size):
                if all(compatible[x, y] for x, y in combinations(xs, 2)):
                    g, _, _ = check_data(xs, scale)
                    maximum_triangle_content = max(maximum_triangle_content, g)
                    checks += 1

    # Arbitrary rational projective parameters, cleared through every pair.
    rng = Random(20260912)
    for size in range(2, 8):
        for _ in range(24):
            base = tuple(sorted(rng.sample(range(2, 36), size)))
            xs, scale = cleared_clique(base)
            g, _, _ = check_data(xs, scale)
            maximum_triangle_content = max(maximum_triangle_content, g)
            checks += 1

    # The exact four-point Pell counterfamily.
    u, v = 1, 0
    pell_checks = 0
    for exponent in range(1, 5):
        u, v = 9 * u + 20 * v, 4 * u + 9 * v
        assert u * u - 5 * v * v == 1
        scale = 24
        xs = (
            90 * u * v + 30 * v * v + 3,
            96 * u * v,
            160 * v * v + 128 * u * v + 16,
        )
        g, _, _ = check_data(xs, scale)
        maximum_triangle_content = max(maximum_triangle_content, g)
        pell_checks += 1
        checks += 1

    print(f"Affine-content identities checked on {checks} integer cliques.")
    print(f"Pell instances checked: {pell_checks}.")
    print(f"Largest primitive triangle content encountered: {maximum_triangle_content}.")
    print("Verified g | 4L^3 and gcd(g,N) | L, with the sharper prime cases.")


if __name__ == "__main__":
    main()
