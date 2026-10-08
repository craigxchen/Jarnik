"""Exact checks for positive Bezoutians and balanced-prime invisibility.

The general arguments are in integer_cotangent_positive_bezoutian_obstruction.md.
All radius checks use the literal Gaussian denominator lcm and all edge norms.
"""
from itertools import combinations
from math import gcd, isqrt, lcm, prod

from check_integer_cotangent_normalization import (
    edge_norm, edge_quotient, primitive_tuple,
)
from check_least_radius_formula import (
    conj, exact_div, gcd_gaussian, lcm_all, norm,
)


def determinant(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    if not n:
        return 1
    previous = sign = 1
    for col in range(n - 1):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[pivot], a[col] = a[col], a[pivot]
            sign = -sign
        value = a[col][col]
        for i in range(col + 1, n):
            for j in range(col + 1, n):
                numerator = value * a[i][j] - a[i][col] * a[col][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][col] = 0
        previous = value
    return sign * a[-1][-1]


def polynomial(roots):
    # Ascending powers: product (z+x).
    coefficients = [1]
    for x in roots:
        new = [0] * (len(coefficients) + 1)
        for j, c in enumerate(coefficients):
            new[j] += x * c
            new[j + 1] += c
        coefficients = new
    return coefficients


def evaluate(coefficients, x):
    value = 0
    for c in reversed(coefficients):
        value = value * x + c
    return value


def bezoutian(roots):
    p = polynomial(roots)
    e, o = p[::2], p[1::2]
    m = len(roots) // 2
    ee, oo = e + [0] * (m + 1 - len(e)), o + [0] * (m + 1 - len(o))
    b = [[0] * m for _ in range(m)]
    for high in range(1, m + 1):
        for low in range(high):
            c = ee[high] * oo[low] - ee[low] * oo[high]
            for h in range(high - low):
                b[high - 1 - h][low + h] += c
    assert b == [list(row) for row in zip(*b)]
    return p, e, o, b


def resultant(f, g):
    m, n = len(f) - 1, len(g) - 1
    ff, gg = list(reversed(f)), list(reversed(g))
    rows = []
    for shift in range(n):
        rows.append([0] * shift + ff + [0] * (n - shift - 1))
    for shift in range(m):
        rows.append([0] * shift + gg + [0] * (m - shift - 1))
    return determinant(rows)


def check_positive(roots, check_resultants=False):
    k, m = len(roots), len(roots) // 2
    p, e, o, b = bezoutian(roots)
    sums = prod(x + y for x, y in combinations(roots, 2))
    assert determinant(b) == sums > 0
    for size in range(1, m + 1):
        assert determinant([row[:size] for row in b[:size]]) > 0
    for s, t in ((0, 1), (-2, 3), (3, -5)):
        kernel = sum(b[i][j] * s**i * t**j
                     for i in range(m) for j in range(m))
        assert (s - t) * kernel == evaluate(e, s) * evaluate(o, t) \
            - evaluate(e, t) * evaluate(o, s)
    for scale in (2, 7):
        bb = bezoutian([scale * x for x in roots])[3]
        assert all(bb[i][j] == b[i][j] * scale**(2*k - 3 - 2*i - 2*j)
                   for i in range(m) for j in range(m))
    if check_resultants:
        negative = [c * (-1)**i for i, c in enumerate(p)]
        o_squared = [0] * (2 * len(o) - 1)
        o_squared[::2] = o
        assert resultant(e, o) == (-1)**(m*(m-1)//2) * sums
        assert resultant(p, o_squared) == sums**2
        assert resultant(p, negative) == 2**k * prod(roots) * sums**2


def valuation(n, p):
    assert n
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def root_minus_one(p, depth):
    a = next(a for a in range(1, p) if (a*a + 1) % p == 0)
    modulus = p
    for _ in range(1, depth):
        a = next(a + c*modulus for c in range(p)
                 if ((a + c*modulus)**2 + 1) % (p*modulus) == 0)
        modulus *= p
    return a


def denominator(y):
    h = (y, 1)
    return exact_div(conj(h), gcd_gaussian(h, conj(h)))


def check_family(s, p, depth):
    assert s >= 3 and prime(p) and p % 4 == 1 and p > s*s + 1
    power = p**depth
    a0 = root_minus_one(p, depth)
    a = next(a0 + c*power for c in range(1, p + 1)
             if all(valuation((a0 + (c+j)*power)**2 + 1, p) == depth
                    for j in range(s)))
    ys = [a + j*power for j in range(s)] + [power + b for b in range(1, s)]
    assert len(set(ys)) == 2*s - 1
    reduced_denominators = [
        abs(x-y) // gcd(abs(x-y), x*y + 1) for x, y in combinations(ys, 2)
    ]
    scale = lcm(*reduced_denominators)
    assert scale % p and min(ys) == power + 1
    assert scale <= lcm(*range(1, s)) * ((p+s)*power)**(s*(s-1))
    xs = [scale * y for y in ys]
    assert all(x >= scale for x in xs)
    edges = xs + [edge_quotient(x, y, scale) for x, y in combinations(xs, 2)]
    all_edge_n = lcm(*(edge_norm(q, scale) for q in edges))
    rows = primitive_tuple(ys, 1)
    assert all_edge_n == norm(rows[0])
    assert all_edge_n == norm(lcm_all([denominator(y) for y in ys]))
    assert valuation(all_edge_n, p) == depth

    # Every p-layer is precisely the s inside rows versus the anchor and
    # s-1 outside rows. This checks every edge, including anchor edges.
    for i, y in enumerate(ys):
        assert valuation(edge_norm(scale*y, scale), p) == (depth if i < s else 0)
    for i, j in combinations(range(len(ys)), 2):
        q = edge_quotient(xs[i], xs[j], scale)
        expected = depth if (i < s) != (j < s) else 0
        assert valuation(edge_norm(q, scale), p) == expected

    check_positive(ys)
    b = bezoutian(xs)[3]
    assert determinant([[entry % p for entry in row] for row in b]) % p
    outside = [denominator(power + j) for j in range(1, s)]
    outside_lcm = lcm_all(outside)
    for j, d in enumerate(outside, start=1):
        assert 2 * norm(d) >= (power + j)**2 + 1
        assert norm(d) % p
    for i, j in combinations(range(s-1), 2):
        assert norm(gcd_gaussian(outside[i], outside[j])) <= (i-j)**2
    assert all_edge_n % (power * norm(outside_lcm)) == 0
    constant = 2**(s-1) * prod((c-b)**2 for b, c in combinations(range(1, s), 2))
    assert constant * all_edge_n >= power**(2*s-1)
    # Exact squared-radius endpoint comparison, with no floating-point angles.
    assert 16 * constant * all_edge_n >= power**(2*s-5) * min(ys)**4


def main():
    positive_checks = 0
    for k in range(2, 15):
        for number, roots in enumerate((
                [1] * k, list(range(1, k+1)),
                [1 + (j*j + 3*j) % 11 for j in range(k)],
                [2**j + 1 for j in range(k)])):
            check_positive(roots, check_resultants=(k <= 9 and number < 2))
            positive_checks += 1

    roots = [4, 4, 4, 1, 2]
    p, _, _, b = bezoutian(roots)
    assert b == [[55808, 4192], [4192, 1058]]
    assert [[c % 17 for c in row] for row in b] == [[14, 10], [10, 4]]
    assert determinant(b) % 17 == 7
    descending = list(reversed(p))
    h = [[descending[2*j-i+1] if 0 <= 2*j-i+1 <= 5 else 0
          for j in range(5)] for i in range(5)]
    assert [determinant([row[:r] for row in h[:r]]) % 17
            for r in range(1, 6)] == [15, 4, 13, 7, 12]

    family_checks = 0
    for s in range(3, 8):
        p = 17 if s == 3 else next(
            p for p in range(s*s+2, 10*s*s) if p % 4 == 1 and prime(p)
        )
        for depth in range(1, 5):
            check_family(s, p, depth)
            family_checks += 1
    print(f"PASS: {positive_checks} positive Bezoutians, both parities, determinant,")
    print("principal minors, kernel, homogeneous scaling, and exact resultants;")
    print(f"{family_checks} actual cliques, all-edge/Gaussian lcm agreement, balanced")
    print("prime layers, local Bezoutian invertibility, and non-endpoint radius bounds.")


if __name__ == "__main__":
    main()
