"""Exact paired-offset divisor quotient and literal all-edge normalizations."""

from itertools import combinations
from math import atan, gcd, lcm


def cotangents(xs: tuple[int, ...], length: int) -> tuple[int, ...]:
    result = []
    for a, b in combinations(xs, 2):
        num = a * b + length * length
        den = a - b
        assert num % den == 0
        result.append(num // den)
    return tuple(result)


def edge_norm(q: int, length: int) -> tuple[int, int]:
    content = gcd(abs(q), length)
    a = q // content
    b = length // content
    eps = 2 if a % 2 and b % 2 else 1
    return (a * a + b * b) // eps, eps


def divisor_count(value: int) -> int:
    return sum(value % d == 0 for d in range(1, value + 1))


def data(xs: tuple[int, ...], length: int) -> tuple[int, int, int, int, int]:
    assert len(xs) >= 2 and all(a < b for a, b in zip(xs, xs[1:]))
    cot = cotangents(xs, length)
    a = xs[0]
    s = a * a + length * length
    offsets = tuple(x - a for x in xs[1:])
    assert all(s % d == 0 for d in offsets)
    cooffsets = tuple(s // d for d in offsets)
    hd = gcd(*offsets)
    hc = gcd(*cooffsets)
    assert s % (hd * hc) == 0
    h = s // (hd * hc)
    assert h == lcm(*offsets) // hd
    assert all((d // hd) * (c // hc) == h for d, c in zip(offsets, cooffsets))
    assert len({d // hd for d in offsets}) == len(offsets)
    assert len(offsets) <= divisor_count(h)
    norms = tuple(edge_norm(q, length)[0] for q in (*xs, *cot))
    n = lcm(*norms)
    ga = gcd(a, length)
    na, eps = edge_norm(a, length)
    assert n % na == 0
    assert h * hd * hc == ga * ga * eps * na
    assert h <= 2 * n * ga * ga // (hd * hc)
    return h, n, hd, hc, ga


def check() -> None:
    xs = (85, 182, 210)
    length = 70
    assert cotangents(xs, length) == (-210, -182, -1540)
    norms = tuple(edge_norm(q, length)[0] for q in (*xs, *cotangents(xs, length)))
    assert norms == (485, 97, 5, 5, 97, 485)
    assert data(xs, length) == (12125, 485, 1, 1, 5)
    assert 6.46 < 2 * 485 ** 0.25 * atan(length / xs[0]) < 6.47

    short = (157, 182, 447)
    assert data(short, 1) == (290, 212298125, 5, 17, 1)
    assert 1.53 < 2 * 212298125 ** 0.25 * atan(1 / short[0]) < 1.54

    for scale in (2, 3, 7, 11):
        for original, ell in ((xs, length), (short, 1)):
            h0, n0, hd0, hc0, ga0 = data(original, ell)
            h1, n1, hd1, hc1, ga1 = data(
                tuple(scale * x for x in original), scale * ell
            )
            assert h1 == h0 and n1 == n0
            assert (hd1, hc1, ga1) == (
                scale * hd0, scale * hc0, scale * ga0
            )


if __name__ == "__main__":
    check()
    print("PASS: complementary divisor quotient, all-edge norms, scale invariance, and fixtures")
