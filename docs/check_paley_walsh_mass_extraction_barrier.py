"""Finite exact checks for the Paley Walsh-extraction mass barrier."""

from fractions import Fraction
from itertools import combinations
from math import comb, prod
from random import Random

from check_general_four_row_quartic_slack_obstruction import paley


def physical_signs(q: int, b: int) -> list[list[int]]:
    h, _ = paley(q)
    m, r = q + 1, b * q
    signs = [[h[x][j % q + 1] for j in range(r)] for x in range(m)]
    for x in range(m):
        signs[x][x] *= -1
    return signs


def determinant_check(q: int, exhaustive: bool) -> None:
    b, tau = 5, 100
    m, r = q + 1, b * q
    h, _ = paley(q)
    s = physical_signs(q, b)
    # These distinct rational weights test the interval-perturbation
    # calculation; existence of actual clustered prime logs is analytic.
    w = [Fraction(tau) + Fraction(j, r * (r + 1)) for j in range(r)]
    assert sum(x - tau for x in w) < 1
    rng = Random(625101)
    subsets = (
        (tuple(i for i in range(m) if mask & (1 << i))
         for mask in range(1, 1 << m))
        if exhaustive else
        (tuple(sorted(rng.sample(range(m), rng.randrange(2, m + 1))))
         for _ in range(1000))
    )
    checked = 0
    for rows in subsets:
        k = len(rows)
        if k < 2:
            continue
        baseline_sum_sq = sum(
            sum(h[x][j % q + 1] for x in rows) ** 2
            for j in range(r)
        )
        assert baseline_sum_sq == b * k * (m - k)
        row_sum_sq = [sum(s[x][j] for x in rows) ** 2 for j in range(r)]
        equal_margin_4 = k * r - sum(row_sum_sq)
        assert equal_margin_4 >= (b - 4) * k * (k - 1)
        weighted_margin_4 = sum(w[j] * (k - row_sum_sq[j]) for j in range(r))
        assert weighted_margin_4 >= Fraction((tau - 1) * k * (k - 1))
        checked += 1

    quartics = 0
    for rows in combinations(range(m), 4):
        t = sum(prod(s[x][j] for x in rows) for j in range(r))
        negative_weights = [
            w[j] for j in range(r)
            if prod(s[x][j] for x in rows) == -1
        ]
        negative_weight = sum(negative_weights)
        weighted_t = sum(w[j] * prod(s[x][j] for x in rows)
                         for j in range(r))
        assert negative_weight == (sum(w) - weighted_t) / 2
        excess = max(0, abs(t) - (4 * b + 8))
        assert excess * excess <= 9 * b * b * q
        # The four good columns with a flip on these rows could absorb
        # at most the four largest negative-parity weights. All other
        # negative mass must be exceptional, exactly as in (5).
        exceptional_min = negative_weight - sum(
            sorted(negative_weights, reverse=True)[:4]
        )
        assert exceptional_min >= (sum(w) - weighted_t) / 2 - 4 * max(w)
        quartics += 1
    assert quartics == comb(m, 4)
    if q == 11:
        # Label eight retained rows by F_2^3 and fit every physical
        # column to its closest oriented Walsh character. The average
        # parallelogram argument must hold even for this optimal fit.
        n = 8
        chars = [
            [1 if bin(x & a).count("1") % 2 == 0 else -1 for x in range(n)]
            for a in range(n)
        ]
        best_error = sum(
            w[j] * min(
                min(sum(s[x][j] != c[x] for x in range(n)),
                    sum(s[x][j] != -c[x] for x in range(n)))
                for c in chars
            ) for j in range(r)
        )
        planes = {
            tuple(sorted((x, x ^ u, x ^ v, x ^ u ^ v)))
            for x in range(n)
            for u in range(1, n)
            for v in range(1, n)
            if u != v
        }
        assert len(planes) == 14
        mean_negative_mass = sum(
            sum(w[j] for j in range(r)
                if prod(s[x][j] for x in rows) == -1)
            for rows in planes
        ) / len(planes)
        assert best_error >= Fraction(n, 4) * mean_negative_mass
    print(f"q={q}: checked {checked} subsets and {quartics} quadruples")


def main() -> None:
    determinant_check(11, True)
    determinant_check(19, False)
    print("PASS: exact determinant margins and quartic exceptional mass")


if __name__ == "__main__":
    main()
