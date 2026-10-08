"""Finite part of the all-{1,2}-degree scalar height obstruction.

For m >= 16 the companion note proves the assertion by Cauchy--Schwarz.
This checks all remaining degree patterns, not sampled vectors.
"""

from fractions import Fraction
from math import comb


def data(ones, twos):
    m = ones + twos
    degrees = [1] * ones + [2] * twos
    total = sum(degrees)
    assert total % 2 == 0
    coefficients = [1]
    for degree in degrees:
        nxt = [0] * (len(coefficients) + degree)
        for i, value in enumerate(coefficients):
            for j in range(degree + 1):
                nxt[i + j] += value
        coefficients = nxt
    middle = total // 2
    rank = coefficients[middle] - coefficients[middle - 1]

    subset_counts = [0] * (total + 1)
    for a in range(ones + 1):
        for b in range(twos + 1):
            subset_counts[a + 2 * b] += comb(ones, a) * comb(twos, b)
    assert sum(subset_counts) == 2**m
    assert subset_counts == subset_counts[::-1]
    assert subset_counts[middle] % 2 == 0
    cuts = subset_counts[middle] // 2
    height = total * 2 ** (m - 3) - sum(
        count * max(0, weight - middle)
        for weight, count in enumerate(subset_counts)
    )
    # Independent sign-sum form of the height in the proof.
    absolute_sum = sum(
        count * abs(2 * weight - total)
        for weight, count in enumerate(subset_counts)
    )
    assert 4 * height == total * 2 ** (m - 1) - absolute_sum
    assert 2 * cuts <= 2 ** (m - 1)
    return rank, cuts, height


EXPECTED = {
    4: (Fraction(2), 0, 4),
    5: (Fraction(8, 3), 2, 3),
    6: (Fraction(18, 5), 0, 6),
    7: (Fraction(29, 6), 4, 3),
    8: (Fraction(232, 35), 0, 8),
    9: (Fraction(8), 6, 3),
    10: (Fraction(209, 21), 6, 4),
    11: (Fraction(49, 4), 8, 3),
    12: (Fraction(4200, 323), 8, 4),
    13: (Fraction(4901, 300), 8, 5),
    14: (Fraction(20250, 1229), 10, 4),
    15: (Fraction(23100, 1163), 10, 5),
}


def main():
    checked = 0
    for m, expected in EXPECTED.items():
        candidates = []
        for ones in range(0, m + 1, 2):
            twos = m - ones
            rank, cuts, height = data(ones, twos)
            denominator = min(cuts, rank - 1)
            assert height >= 2 * denominator, (ones, twos)
            if cuts and m >= 5:
                assert height > 2 * denominator, (ones, twos)
            if denominator:
                candidates.append((Fraction(height, denominator), ones, twos))
            checked += 1
        assert min(candidates) == expected, (m, min(candidates), expected)
    assert [data(a, 4 - a) for a in (4, 2, 0)] == [
        (2, 3, 2), (2, 2, 2), (3, 3, 4)
    ]
    print(f"Passed: all {checked} degree patterns for 4 <= m <= 15.")
    print("Passed: exact ratio table, four-row equality, sign-height identity.")
    print("All m >= 16 are covered by the inequality in the companion note.")


if __name__ == "__main__":
    main()
