#!/usr/bin/env python3
"""Exact finite checks for Gaussian content and common-rotation claims."""

from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd, isqrt

Gaussian = tuple[int, int]


def add(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] + w[0], z[1] + w[1]


def sub(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] - w[0], z[1] - w[1]


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def norm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def div_exact(z: Gaussian, w: Gaussian) -> Gaussian:
    n = norm(w)
    a, b = mul(z, conj(w))
    assert a % n == 0 and b % n == 0, (z, w)
    return a // n, b // n


def divmod_gaussian(z: Gaussian, w: Gaussian) -> tuple[Gaussian, Gaussian]:
    n = norm(w)
    a, b = mul(z, conj(w))
    # Nearest integer, with half ties rounded upward; either nearest choice
    # gives the usual strict Gaussian Euclidean remainder bound.
    q = ((2 * a + n) // (2 * n), (2 * b + n) // (2 * n))
    return q, sub(z, mul(q, w))


def gcd_gaussian(z: Gaussian, w: Gaussian) -> Gaussian:
    while w != (0, 0):
        _, r = divmod_gaussian(z, w)
        z, w = w, r
    return z


def gcd_many_gaussian(values: tuple[Gaussian, ...]) -> Gaussian:
    return reduce(gcd_gaussian, values)


def integer_content(rows: tuple[Gaussian, ...]) -> int:
    ns = [norm(z) for z in rows]
    ds = [mul(conj(rows[i]), rows[j])[1] for i in range(len(rows)) for j in range(i + 1, len(rows))]
    return reduce(gcd, ns + ds)


def signature(rows: tuple[Gaussian, ...]) -> tuple[Gaussian, ...]:
    return tuple(
        mul(conj(rows[i]), rows[j])
        for i in range(len(rows))
        for j in range(i, len(rows))
    )


def candidates_of_norm(n: int) -> tuple[Gaussian, ...]:
    return tuple(
        (a, b)
        for a in range(-isqrt(n), isqrt(n) + 1)
        for b in range(-isqrt(n), isqrt(n) + 1)
        if a * a + b * b == n
    )


def gaussian_primitive(z: Gaussian) -> bool:
    # Equivalent to gcd(z, conjugate(z)) being a unit.
    return gcd(abs(z[0]), abs(z[1])) == 1 and (z[0] - z[1]) % 2 != 0


def all_fixed_gram_realizations(rows: tuple[Gaussian, ...]) -> set[tuple[Gaussian, ...]]:
    h = gcd_many_gaussian(rows)
    d = norm(h)
    qs = tuple(div_exact(z, h) for z in rows)
    return {
        tuple(mul(hp, q) for q in qs)
        for hp in candidates_of_norm(d)
    }


def check_exhaustive_small_rows() -> None:
    points = tuple((x, y) for x in range(-2, 3) for y in range(-2, 3) if (x, y) != (0, 0))
    groups: dict[tuple[Gaussian, ...], set[tuple[Gaussian, ...]]] = {}
    tested = 0
    for m in (2, 3):
        for rows in product(points, repeat=m):
            rows = tuple(rows)
            d = norm(gcd_many_gaussian(rows))
            d0 = integer_content(rows)
            assert d == d0, (rows, d, d0)
            groups.setdefault(signature(rows), set()).add(rows)
            tested += 1

    # Every possible row of each fixed Gram signature has norm at most 8,
    # so the square [-2,2]^2 is exhaustive for these classes.
    for sig, observed in groups.items():
        exemplar = next(iter(observed))
        classified = all_fixed_gram_realizations(exemplar)
        assert classified == observed, (sig, len(classified), len(observed))
    print(f"small rows: content identity and rotation classes checked ({tested} tuples, {len(groups)} classes)")


def check_finite_primitive_fixture() -> None:
    rows = ((1_469_178, 1), (31_686, 1), (153_548, 1))
    h = (10, 3)
    d = norm(h)
    assert d == 109
    qs = tuple(div_exact(z, h) for z in rows)
    assert norm(gcd_many_gaussian(rows)) == d
    assert integer_content(rows) == d
    assert all(gaussian_primitive(z) for z in rows)
    assert gcd_many_gaussian(qs) in ((1, 0), (-1, 0), (0, 1), (0, -1))

    b = 1
    # A right-half-circle cap with |Im|<=B has chord diameter at most 2B.
    assert min(norm(z) for z in rows) > 4 * d * b * b
    candidates = all_fixed_gram_realizations(rows)
    assert len(candidates) == 8  # r_2(109)
    narrow_positive = {
        candidate
        for candidate in candidates
        if all(z[0] > 0 and abs(z[1]) <= b for z in candidate)
    }
    assert narrow_positive == {rows}


def check_underfull_profile_bound() -> None:
    common = (1, 1)       # norm 2
    pair12 = (3, 2)       # norm 13
    pair23 = (4, 1)       # norm 17
    pair13 = (5, 2)       # norm 29
    private = (7, 0)  # common row factor K_i, norm 49
    core_norms = (norm(common), norm(pair12), norm(pair23), norm(pair13))
    assert reduce(gcd, core_norms) == 1
    rows = (
        mul(private, mul(common, mul(pair12, pair13))),
        mul(private, mul(common, mul(pair12, pair23))),
        mul(private, mul(common, mul(pair13, pair23))),
    )
    d = norm(gcd_many_gaussian(rows))
    g = norm(common)
    k_i = norm(private)
    assert (d, g, k_i) == (98, 2, 49)
    assert d % g == 0 and (k_i**3) % (d // g) == 0


def check_rational_square_obstruction() -> None:
    # Q=I, v1=(3,4), v2=(-7/5,24/5) have norms 25 and
    # dot/area (15,20), but no integral realization is primitive in both rows.
    v1 = (Fraction(3), Fraction(4))
    v2 = (Fraction(-7, 5), Fraction(24, 5))
    dot = v1[0] * v2[0] + v1[1] * v2[1]
    delta = v1[0] * v2[1] - v1[1] * v2[0]
    assert dot == 15 and delta == 20
    assert sum(x * x for x in v1) == 25
    assert sum(x * x for x in v2) == 25

    z = (15, 20)
    realizations = []
    for g in candidates_of_norm(25):
        try:
            p2 = div_exact(z, g)
        except AssertionError:
            continue
        p1 = conj(g)
        assert norm(p1) == norm(p2) == 25
        assert mul(conj(p1), p2) == z
        realizations.append((p1, p2))
    assert realizations
    assert not any(gaussian_primitive(p1) and gaussian_primitive(p2) for p1, p2 in realizations)


def main() -> None:
    check_exhaustive_small_rows()
    check_finite_primitive_fixture()
    check_underfull_profile_bound()
    check_rational_square_obstruction()
    print("Gaussian content/common-rotation checks passed")


if __name__ == "__main__":
    main()
