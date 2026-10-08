"""Exact cut and Gaussian checks for half-symmetric star isolation."""

from fractions import Fraction
from itertools import combinations
from math import isqrt


def mul(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def conj(z):
    return z[0], -z[1]


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def split_primes(n):
    out = []
    p = 5
    while len(out) < n:
        if p % 4 == 1 and all(p % d for d in range(2, isqrt(p) + 1)):
            out.append(p)
        p += 2
    return out


def gaussian_prime(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b and b * b == b2:
            return a, b
    raise AssertionError(p)


def check_cut_counts():
    checked = 0
    for k in range(2, 7):
        a = set(range(k))
        b = set(range(k, 2 * k))
        alpha = Fraction(k, 2 * (k - 1))
        sharp = Fraction((k * k) // 4, k * (k - 1) // 2)
        for t in range(1, k):
            for u in combinations(a, t):
                for v in combinations(b, t):
                    s = set(u + v)
                    assert len(s & a) == len(s & b) == t
                    assert not (s <= a or a <= s)
                    cross = sum((i in s) != (j in s) for i in a for j in b)
                    within = sum((i in s) != (j in s) for i, j in combinations(a, 2))
                    assert cross == 2 * t * (k - t)
                    assert within == t * (k - t)
                    assert Fraction(cross, k * k) <= Fraction(1, 2)
                    assert Fraction(within, k * (k - 1) // 2) <= sharp <= alpha
                    checked += 1
        assert sharp == (Fraction(2, 3) if k == 3 else sharp)
        assert sharp == alpha if k % 2 == 0 else sharp < alpha
    return checked


def check_actual_gaussian_tuple():
    k = 5
    a = list(range(k))
    b = list(range(k, 2 * k))
    regular = [set((i, j)) for i in a for j in b]  # t=1; not half-balanced.
    cuts = regular + [set(a)]
    primes = split_primes(len(cuts) + 1)
    bases = [gaussian_prime(p) for p in primes[:len(cuts)]]
    exponents = [1 + (j % 3) for j in range(len(cuts))]
    star = power(bases[-1], exponents[-1])
    d = gaussian_prime(primes[-1])  # A complete common split multiplier.
    epsilon = (0, 1)
    rows = []
    for i in range(2 * k):
        z = mul(epsilon, d)
        for cut, pi, exponent in zip(cuts, bases, exponents):
            z = mul(z, power(pi if i in cut else conj(pi), exponent))
        rows.append(z)
    n = norm(rows[0])
    assert len(set(rows)) == 2 * k
    assert all(norm(z) == n for z in rows)
    prod_a = prod_b = (1, 0)
    for i in a:
        prod_a = mul(prod_a, rows[i])
    for j in b:
        prod_b = mul(prod_b, rows[j])
    assert mul(prod_a, power(conj(star), k)) == mul(prod_b, power(star, k))

    # The literal common unit and multiplier survive each pair factorization.
    pairs = 0
    for i, j in combinations(range(2 * k), 2):
        pair = (1, 0)
        common = mul(epsilon, d)
        for cut, pi, exponent in zip(cuts, bases, exponents):
            if (i in cut) == (j in cut):
                common = mul(common, power(pi if i in cut else conj(pi), exponent))
            else:
                pair = mul(pair, power(pi if i in cut else conj(pi), exponent))
        assert mul(common, pair) == rows[i]
        assert mul(common, conj(pair)) == rows[j]
        assert norm(pair) * norm(common) == n
        dx, dy = rows[i][0] - rows[j][0], rows[i][1] - rows[j][1]
        assert dx * dx + dy * dy == 4 * pair[1] * pair[1] * norm(common)
        pairs += 1
    assert pairs == 45
    return len(cuts), pairs


if __name__ == "__main__":
    cuts, (source_columns, pairs) = check_cut_counts(), check_actual_gaussian_tuple()
    print(f"PASS: {cuts} equal-half cut incidences through k=6; a literal common-unit, "
          f"common-multiplier ten-row tuple with {source_columns} source columns, "
          f"isolated star phase, and {pairs} exact pair chord identities.")
