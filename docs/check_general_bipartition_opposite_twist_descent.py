"""Exact checks for arbitrary-bipartition opposite-twist descent."""

from functools import reduce
from fractions import Fraction
from itertools import product
from math import comb, prod
from random import Random

from check_ordered_affine_circuit_conductor_index import (
    gconj,
    gdiv_exact,
    ggcd,
    gmul,
    gnorm,
    gpower,
)


UNITS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def cuts(m):
    """Unordered nontrivial cuts, represented with the last row outside."""
    universe = set(range(m))
    for mask in range(1, 1 << (m - 1)):
        side = {j for j in range(m - 1) if mask >> j & 1}
        yield side, universe - side


def oriented_gap(t, side, other):
    lo_s, hi_s = min(t[j] for j in side), max(t[j] for j in side)
    lo_t, hi_t = min(t[j] for j in other), max(t[j] for j in other)
    if lo_s > hi_t:
        return lo_s - hi_t, 1
    if lo_t > hi_s:
        return lo_t - hi_s, -1
    return 0, 0


def valuation(z, pi):
    answer = 0
    while True:
        try:
            z = gdiv_exact(z, pi)
        except AssertionError:
            return answer
        answer += 1


def descend(points, prime_data, side, other):
    alpha = (1, 0)
    k = 1
    for pi, e, t in prime_data:
        h, orientation = oriented_gap(t, side, other)
        factor = pi if orientation == 1 else gconj(pi)
        alpha = gmul(alpha, gpower(factor, h))
        k *= gnorm(pi) ** h
    descended = [
        gdiv_exact(z, alpha if j in side else gconj(alpha))
        for j, z in enumerate(points)
    ]
    return alpha, k, descended


def descend_allocations(t, e, side, other):
    h, orientation = oriented_gap(t, side, other)
    if orientation == 1:
        out = [v - h if j in side else v for j, v in enumerate(t)]
    elif orientation == -1:
        out = [v if j in side else v - h for j, v in enumerate(t)]
    else:
        out = list(t)
    return h, out, e - h


def main():
    allocation_tables = 0
    for m in range(2, 7):
        for e in range(1, 6):
            for t in product(range(e + 1), repeat=m):
                if min(t) != 0 or max(t) != e:
                    continue
                assert sum(oriented_gap(t, s, u)[0] for s, u in cuts(m)) == e
                if m <= 5:
                    all_cuts = list(cuts(m))
                    old = [oriented_gap(t, s, u) for s, u in all_cuts]
                    for chosen, (s, u) in enumerate(all_cuts):
                        if old[chosen][0] == 0:
                            continue
                        h, after, enew = descend_allocations(t, e, s, u)
                        assert h == old[chosen][0]
                        assert min(after) == 0 and max(after) == enew
                        new = [
                            oriented_gap(after, v, w) for v, w in all_cuts
                        ]
                        assert new[chosen][0] == 0
                        assert all(
                            new[j] == old[j]
                            for j in range(len(all_cuts))
                            if j != chosen
                        )
                allocation_tables += 1

    for m in range(3, 21):
        rho = min(
            Fraction(comb(s, 3) + comb(m - s, 3), comb(m, 3))
            for s in range(1, m)
        )
        k = m // 2
        expected = (
            Fraction(k - 2, 2 * (2 * k - 1))
            if m == 2 * k
            else Fraction(k - 1, 2 * (2 * k + 1))
        )
        assert rho == expected < Fraction(1, 4)

    rng = Random(923611)
    pis = ((2, 1), (3, 2), (4, 1))
    random_tuples = 0
    strict_descents = 0
    for _ in range(500):
        m = rng.randrange(3, 8)
        points = [rng.choice(UNITS) for _ in range(m)]
        prime_data = []
        n = 1
        for pi in pis:
            e0 = rng.randrange(1, 6)
            raw = [rng.randrange(e0 + 1) for _ in range(m)]
            lo, hi = min(raw), max(raw)
            if lo == hi:
                continue
            t = [v - lo for v in raw]
            e = hi - lo
            prime_data.append((pi, e, t))
            n *= gnorm(pi) ** e
            for j in range(m):
                points[j] = gmul(
                    points[j],
                    gmul(gpower(pi, t[j]), gpower(gconj(pi), e - t[j])),
                )
        assert all(gnorm(z) == n for z in points)
        all_cuts = list(cuts(m))
        gap_data = []
        for side, other in all_cuts:
            k = prod(
                gnorm(pi) ** oriented_gap(t, side, other)[0]
                for pi, _, t in prime_data
            )
            gap_data.append((k, side, other))
        assert prod(k for k, _, _ in gap_data) == n
        k_expected, side, other = max(gap_data, key=lambda row: row[0])
        alpha, k, out = descend(points, prime_data, side, other)
        assert k == k_expected == gnorm(alpha)
        assert k > 1
        assert all(gnorm(z) == n // k for z in out)
        assert reduce(ggcd, out) in UNITS
        assert all(
            oriented_gap([valuation(z, pi) for z in out], side, other)[0] == 0
            for pi, _, _ in prime_data
        )
        gcd_s = reduce(ggcd, (points[j] for j in side))
        gcd_t = reduce(ggcd, (points[j] for j in other))
        assert gnorm(gcd_s) % k == 0 and gnorm(gcd_t) % k == 0
        strict_descents += k > 1
        random_tuples += 1

    # The second-stage exact collision in (22).
    stage_one = [(-14, 5), (11, 10), (10, 11), (5, 14), (-14, -5)]
    side, other = {0, 2, 3}, {1, 4}
    pi13, pi17 = (3, 2), (4, 1)
    t13 = [valuation(z, pi13) for z in stage_one]
    t17 = [valuation(z, pi17) for z in stage_one]
    assert t13 == [0, 0, 1, 0, 1]
    assert t17 == [1, 0, 1, 1, 0]
    assert oriented_gap(t13, side, other)[0] == 0
    alpha, k, stage_two = descend(
        stage_one, [(pi13, 1, t13), (pi17, 1, t17)], side, other
    )
    assert alpha == (4, 1) and k == 17
    assert all(gnorm(z) == 13 for z in stage_two)
    assert stage_two[1] == stage_two[3] == (2, 3)
    # Collision distance identity (9): |z_3-z_1|^2=4 b^2 M.
    dx = stage_one[3][0] - stage_one[1][0]
    dy = stage_one[3][1] - stage_one[1][1]
    assert dx * dx + dy * dy == 4 * alpha[1] ** 2 * 13

    # Across both descents, all original multipliers have norm Q=85 and
    # are congruent modulo 2.  The colliding rows satisfy (14)--(15).
    alpha1, alpha2 = (-1, 2), (4, 1)
    beta1 = gmul(gconj(alpha1), gconj(alpha2))
    beta3 = gmul(gconj(alpha1), alpha2)
    assert gnorm(beta1) == gnorm(beta3) == 85
    assert (beta1[0] - beta3[0]) % 2 == 0
    assert (beta1[1] - beta3[1]) % 2 == 0
    w = stage_two[1]
    assert gmul(beta1, w) == (9, -32)
    assert gmul(beta3, w) == (23, -24)

    print(
        f"PASS: {allocation_tables} allocation tables, {random_tuples} "
        "random all-cut products/descents, and the two-stage collision fixture"
    )


if __name__ == "__main__":
    main()
