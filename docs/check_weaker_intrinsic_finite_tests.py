"""Exact audits for weaker_intrinsic_finite_tests.md; no endpoint search."""

from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import comb


def mean(values):
    k = len(values) - 1
    return sum(Fraction(comb(k, r)) * x for r, x in enumerate(values)) / 2**k


def check_binomial_tests():
    expected = {
        6: Fraction(3, 4),
        8: Fraction(2, 5),
        10: Fraction(5, 28),
        12: Fraction(1, 14),
        14: Fraction(7, 264),
        16: Fraction(4, 429),
    }
    for k, threshold in expected.items():
        assert threshold == Fraction(k * (k - 1), 2 * comb(k, k // 2))
        for delta in (threshold, Fraction(1, 2), Fraction(1)):
            values = [
                Fraction(k, 4)
                if r in (0, k)
                else Fraction((r - k // 2) ** 2) + delta * (r == k // 2)
                for r in range(k + 1)
            ]
            gap = (delta * comb(k, k // 2) - Fraction(k * (k - 1), 2)) / 2**k
            assert mean(values) - Fraction(k, 4) == gap
    for a, b in ((8, 3), (6, 3), (2, 4), (0, 5), (8, Fraction(5, 2))):
        values = [2, a + 1, b + 1, 1, 1, 1, b + 1, a + 1, 2]
        assert mean(values) - 2 == Fraction(8 * a + 28 * b - 127, 128)
    assert 42 * Fraction(17, 2) / Fraction(7, 256) == 13056
    assert 42 * 6 / Fraction(5, 128) == Fraction(32256, 5)


def check_source_gcd_dictionary():
    cases = 0
    for allocations in combinations_with_replacement(range(5), 8):
        if allocations[0] != 0 or allocations[-1] == 0:
            continue
        e = allocations[-1]
        v = [0] * 5
        for ell in range(1, e + 1):
            r = sum(a >= ell for a in allocations)
            v[min(r, 8 - r)] += 1
        g = []
        for deleted_size in range(1, 4):
            total = 0
            for deleted in combinations(range(8), deleted_size):
                remaining = [a for i, a in enumerate(allocations) if i not in deleted]
                # Contributions of the two Gaussian primes over the rational prime.
                total += min(remaining) + e - max(remaining)
            g.append(total)
        assert g[0] == v[1]
        assert g[1] == 7 * v[1] + v[2]
        assert g[2] == 21 * v[1] + 6 * v[2] + v[3]
        assert e == sum(v)
        # Both exact half-bonus source-gcd reformulations.
        half_excess_twice = 14 * v[1] + 4 * v[2] - 2 * v[3] - 3 * v[4]
        assert half_excess_twice == 16 * v[1] + 6 * v[2] - 2 * e - v[4]
        assert half_excess_twice == g[1] + g[2] - 3 * e - 11 * g[0]
        # The reduced-singleton target in the two stated forms.
        excess = 5 * v[1] + 2 * v[2] - v[3] - v[4]
        assert excess == 6 * v[1] + 3 * v[2] - e
        assert excess == 3 * g[1] - e - 15 * g[0]
        cases += 1
    assert cases == 329
    return cases


if __name__ == "__main__":
    check_binomial_tests()
    count = check_source_gcd_dictionary()
    print(f"Passed exact binomial tests and {count} prime-power source-gcd audits.")
