"""Exact checks for the product of source-anchor-dependent diagonal radii."""

from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, prod
from random import Random


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def det(z, w):
    return z[0] * w[1] - z[1] * w[0]


def lcm(x, y):
    return x // gcd(x, y) * y


def ggcd(z, w):
    while w != (0, 0):
        numerator = mul(z, conj(w))
        denominator = norm(w)
        q = tuple((2 * n + denominator) // (2 * denominator)
                  for n in numerator)
        p = mul(q, w)
        z, w = w, (z[0] - p[0], z[1] - p[1])
    return z


def edge_norm(z):
    g = gcd(abs(z[0]), abs(z[1]))
    x, y = z[0] // g, z[1] // g
    return (x * x + y * y) // (2 if x % 2 and y % 2 else 1)


def radius(rows):
    out = 1
    for z, w in combinations(rows, 2):
        out = lcm(out, edge_norm(mul(z, conj(w))))
    return out


def away(x, support):
    while gcd(x, support) != 1:
        x //= gcd(x, support)
    return x


def vp(n, p):
    assert n > 0
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


@lru_cache(None)
def split_prime(p):
    assert p % 4 == 1
    for x in range(1, isqrt(p) + 1):
        y = isqrt(p - x * x)
        if x * x + y * y == p:
            return x, y
    raise AssertionError(p)


def vg(z, pi):
    assert z != (0, 0)
    p = norm(pi)
    out = 0
    while True:
        numerator = mul(z, conj(pi))
        if numerator[0] % p or numerator[1] % p:
            return out
        z = numerator[0] // p, numerator[1] // p
        out += 1


def check(rows, a, d):
    assert rows[0] == (1, 0)
    assert all(gcd(abs(x), abs(y)) == 1 for x, y in rows)
    assert all(det(z, w) for z, w in combinations(rows, 2))
    assert gcd(a, d) == 1 and a > 0 and d > 0
    m = len(rows)
    source_n = radius(rows)
    b = {(i, k): abs(det(rows[i], rows[k])) // norm(ggcd(rows[i], rows[k]))
         for i, k in combinations(range(m), 2)}
    big_b = prod(b.values())
    targets, target_n, edges = [], [], {}
    u, v = a + d, a - d
    for j in range(m):
        transformed = []
        for i in range(m):
            w = mul(rows[i], conj(rows[j]))
            y = a * w[0], d * w[1]
            # Exact rational-circle Mobius numerator after common clearing.
            assert (u * w[0] + v * w[0], u * w[1] - v * w[1]) == \
                (2 * y[0], 2 * y[1])
            transformed.append(y)
            if i < j:
                edges[i, j] = edge_norm(y)
        targets.append(transformed)
        target_n.append(radius(transformed))
    numerator = prod(n * n for n in edges.values())
    denominator = prod(target_n)
    assert numerator % denominator == 0
    quotient = numerator // denominator
    assert all(n % 2 for n in [source_n, quotient] + list(edges.values()) + target_n)
    out_q = away(quotient, 2 * source_n)
    out_b = away(big_b, 2 * source_n)
    assert out_b ** (m - 2) % out_q == 0
    assert away(quotient, 2 * source_n * big_b) == 1

    for p in (5, 13, 17, 29, 37, 41, 53):
        pi = split_prime(p)
        t = [[vg(targets[j][i], pi) - vg(targets[j][i], conj(pi))
              for j in range(m)] for i in range(m)]
        assert all(t[i][j] == -t[j][i] for i in range(m) for j in range(m))
        losses = []
        for j in range(m):
            col = [t[i][j] for i in range(m)]
            assert vp(target_n[j], p) == max(col) - min(col)
            losses.append(sum(abs(x) for x in col) - max(col) + min(col))
            if source_n % p:
                for i, k in combinations([i for i in range(m) if i != j], 2):
                    if t[i][j] * t[k][j] > 0:
                        assert min(abs(t[i][j]), abs(t[k][j])) <= vp(b[i, k], p)
                assert losses[-1] <= sum(vp(value, p) for (i, k), value in b.items()
                                        if i != j and k != j)
        assert sum(losses) == vp(quotient, p)
        assert all(vp(edges[i, j], p) == abs(t[i][j]) for i, j in edges)
        if source_n % p and (u % p == 0 or v % p == 0):
            assert all(x == 0 for row in t for x in row)
    return source_n, big_b, list(edges.values()), target_n, quotient


def main():
    first = check([(1, 0), (2, 1), (3, 1)], 2, 1)
    assert first[2:] == ([17, 37, 197], [629, 3349, 7289], 1)
    second = check([(1, 0), (4, 1), (9, 1)], 1, 2)
    assert second == (697, 5, [5, 85, 1469], [85, 7345, 124865], 5)
    rng = Random(572960)
    matrices = ((1, 1), (2, 1), (1, 2), (3, 2), (4, 1),
                (7, 2), (5, 1), (1, 5), (13, 2), (2, 13))
    fixtures = [
        [(1, 0), (4, 1), (9, 1)],
        [(1, 0), (1, 1), (2, 1), (1, 2)],
        [(1, 0), (4, 1), (9, 1), (14, 1), (19, 1)],
    ]
    for _ in range(50):
        rows = [(1, 0)]
        size = rng.randrange(3, 7)
        while len(rows) < size:
            row = (rng.randrange(-10, 11), rng.randrange(1, 11))
            if gcd(abs(row[0]), row[1]) == 1 and all(det(row, h) for h in rows):
                rows.append(row)
        fixtures.append(rows)
    count = 2
    nontrivial = 0
    for rows in fixtures:
        for a, d in matrices:
            n, b, edges, targets, quotient = check(rows, a, d)
            nontrivial += away(quotient, 2 * n) > 1
            count += 1
    print(f"PASS: {count} exact anchor-indexed radius-product cases; "
          f"{nontrivial} nontrivial residue corrections away from the source radius; "
          "seven split-prime local audits per case; parity retained.")


if __name__ == "__main__":
    main()
