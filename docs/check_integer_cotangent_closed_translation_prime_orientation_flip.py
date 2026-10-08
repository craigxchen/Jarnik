"""Exact fixture for closed cotangent translation and cut-prime orientation."""

from itertools import combinations
from math import atan, gcd, lcm


def cotangents(xs: tuple[int, ...], length: int) -> dict[tuple[int, int], int]:
    result = {}
    for i, j in combinations(range(len(xs)), 2):
        num = xs[i] * xs[j] + length * length
        den = xs[i] - xs[j]
        assert num % den == 0
        result[i, j] = num // den
    return result


def edge_norm(q: int, length: int) -> int:
    common = gcd(abs(q), length)
    a = q // common
    b = length // common
    eps = 2 if a % 2 and b % 2 else 1
    return (a * a + b * b) // eps


def radius_squared(xs: tuple[int, ...], length: int) -> int:
    qs = cotangents(xs, length)
    return lcm(*(edge_norm(q, length) for q in (*xs, *qs.values())))


def valuation(value: int, prime: int) -> int:
    count = 0
    while value % prime == 0:
        count += 1
        value //= prime
    return count


def check_fixture() -> None:
    length = 1
    source = (31, 32, 57)
    shift = 250
    target = tuple(x + shift for x in source)
    assert source == (31, 32, 57)
    assert target == (281, 282, 307)
    assert tuple(cotangents(source, length).values()) == (-993, -68, -73)
    assert tuple(cotangents(target, length).values()) == (-79243, -3318, -3463)

    deltas = [abs(a - b) for a, b in combinations(source, 2)]
    assert deltas == [1, 26, 25]
    assert lcm(*deltas) == 650
    assert shift % 13 == 3  # -2*5 modulo 13
    assert shift % 25 == 0
    assert shift % 2 == 0
    for a, b in combinations(source, 2):
        assert shift * (a + b + shift) % (a - b) == 0

    assert tuple(x % 13 for x in source) == (5, 6, 5)
    assert tuple(x % 13 for x in target) == (8, 9, 8)
    assert {r for r in range(13) if (r * r + 1) % 13 == 0} == {5, 8}
    assert tuple(x % 5 for x in source) == (1, 2, 2)
    assert tuple(x % 5 for x in target) == (1, 2, 2)
    assert {r for r in range(5) if (r * r + 1) % 5 == 0} == {2, 3}

    n0 = radius_squared(source, length)
    n1 = radius_squared(target, length)
    assert (n0, n1) == (2465125, 455260346125)
    assert valuation(n0, 13) == valuation(n1, 13) == 1
    assert valuation(n0, 5) == valuation(n1, 5) == 3
    c0 = 2 * n0 ** 0.25 * atan(length / min(source))
    c1 = 2 * n1 ** 0.25 * atan(length / min(target))
    assert 2.55 < c0 < 2.56
    assert 5.84 < c1 < 5.85


def check_two_independent_flips() -> None:
    source = (157, 182, 447)
    assert tuple(cotangents(source, 1).values()) == (-1143, -242, -307)
    assert lcm(*(abs(a - b) for a, b in combinations(source, 2))) == 76850
    cases = {
        (False, False): (0, 212298125, (12, 8, 12), (51, 23, 23)),
        (False, True): (
            65250, 102866995511848031667125, (12, 8, 12), (5, 30, 30)
        ),
        (True, False): (
            29150, 842320937321019957125, (17, 13, 17), (51, 23, 23)
        ),
        (True, True): (
            17550, 41549115737611768325, (17, 13, 17), (5, 30, 30)
        ),
    }
    assert {r for r in range(29) if (r * r + 1) % 29 == 0} == {12, 17}
    assert {r for r in range(53) if (r * r + 1) % 53 == 0} == {23, 30}
    for (flip29, flip53), (shift, expected_n, residues29, residues53) in cases.items():
        assert shift % 50 == 0
        assert shift % 29 == (5 if flip29 else 0)
        assert shift % 53 == (7 if flip53 else 0)
        xs = tuple(x + shift for x in source)
        assert tuple(x % 29 for x in xs) == residues29
        assert tuple(x % 53 for x in xs) == residues53
        for a, b in combinations(source, 2):
            assert shift * (a + b + shift) % (a - b) == 0
        n = radius_squared(xs, 1)
        assert n == expected_n
        assert valuation(n, 29) == valuation(n, 53) == 1
        c = 2 * n ** 0.25 * atan(1 / min(xs))
        if not flip29 and not flip53:
            assert 1.53 < c < 1.54
        else:
            assert c > 9


if __name__ == "__main__":
    check_fixture()
    check_two_independent_flips()
    print("PASS: exact closed translation, independent cut-prime orientation, and all-edge radii")
