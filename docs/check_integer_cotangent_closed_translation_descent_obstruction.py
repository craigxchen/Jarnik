"""Exact CRT and all-edge certificate for the translation descent obstruction."""

from itertools import combinations
from math import lcm


SOURCE = (157, 182, 447)
CLASSES = (0, 17550, 29150, 31436, 43036, 60586, 65250, 72186)
MODULUS = 76850


def closed(shift: int) -> bool:
    return all(
        shift * (a + b + shift) % (a - b) == 0
        for a, b in combinations(SOURCE, 2)
    )


def cotangents(xs: tuple[int, ...]) -> tuple[int, ...]:
    result = []
    for a, b in combinations(xs, 2):
        num = a * b + 1
        den = a - b
        assert num % den == 0
        result.append(num // den)
    return tuple(result)


def edge_norm(q: int) -> int:
    eps = 2 if q % 2 else 1
    return (q * q + 1) // eps


def all_edge_radius_squared(xs: tuple[int, ...]) -> int:
    return lcm(*(edge_norm(q) for q in (*xs, *cotangents(xs))))


def check() -> None:
    assert cotangents(SOURCE) == (-1143, -242, -307)
    assert lcm(*(abs(a - b) for a, b in combinations(SOURCE, 2))) == MODULUS
    assert all_edge_radius_squared(SOURCE) == 212298125
    # Rational certificates for 1.53<C_*<1.54, using
    # u-u^3/3 < atan(u) < u and raising positive quantities to power four.
    a, n = min(SOURCE), 212298125
    assert 16*n*(3*a*a-1)**4*100**4 > 153**4*(3*a**3)**4
    assert 16*n*100**4 < 154**4*a**4

    assert tuple(t for t in range(MODULUS) if closed(t)) == CLASSES
    assert all(t % 2 == 0 for t in CLASSES)
    assert {t % 25 for t in CLASSES} == {0, 11}
    assert {t % 29 for t in CLASSES} == {0, 5}
    assert {t % 53 for t in CLASSES} == {0, 7}
    assert min(min(t, MODULUS - t) for t in CLASSES if t) == 4664

    pair_lower = (4507 * 4482 + 1) // 25
    assert pair_lower == 808015
    assert pair_lower * pair_lower // 2 > 212298125
    for t in (-MODULUS - 4664, -4664, 17550, MODULUS + 17550):
        assert closed(t)
        xs = tuple(x + t for x in SOURCE)
        q12 = cotangents(xs)[0]
        assert abs(q12) >= pair_lower
        assert edge_norm(q12) > 212298125
        assert all_edge_radius_squared(xs) > 212298125


if __name__ == "__main__":
    check()
    print("PASS: eight exact CRT classes and all nonzero closed translations exceed source radius")
