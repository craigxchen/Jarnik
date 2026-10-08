"""Exact dimension, weighted-cut, and height checks for degree variation."""

from fractions import Fraction
from itertools import product
from math import comb


def data(degrees):
    total = sum(degrees)
    assert total % 2 == 0
    coefficients = [1]
    for d in degrees:
        nxt = [0] * (len(coefficients) + d)
        for i, c in enumerate(coefficients):
            for j in range(d + 1):
                nxt[i + j] += c
        coefficients = nxt
    r = coefficients[total // 2] - coefficients[total // 2 - 1]
    cut_weights = [sum(d for d, bit in zip(degrees, bits) if bit)
                  for bits in product((0, 1), repeat=len(degrees))]
    s = sum(2 * x == total for x in cut_weights) // 2
    height = total * 2 ** (len(degrees) - 3)
    height -= sum(max(0, x - total // 2) for x in cut_weights)
    return r, s, height


expected = {
    (4, 1): (2, 2, 2, Fraction(2)),
    (4, 2): (3, 3, 4, Fraction(2)),
    (4, 3): (4, 3, 6, Fraction(2)),
    (4, 4): (5, 3, 8, Fraction(8, 3)),
    (6, 1): (5, 5, 18, Fraction(9, 2)),
    (6, 2): (15, 10, 36, Fraction(18, 5)),
    (6, 3): (34, 10, 54, Fraction(27, 5)),
    (6, 4): (65, 10, 72, Fraction(36, 5)),
    (8, 1): (14, 14, 116, Fraction(116, 13)),
    (8, 2): (91, 35, 232, Fraction(232, 35)),
    (8, 3): (364, 35, 348, Fraction(348, 35)),
}
for (m, d), expected_row in expected.items():
    r, s, height = data([d] * m)
    rank = comb(m, m // 2) // (m // 2 + 1) if d == 1 else s
    denominator = rank if r > rank else r - 1
    got = (r, rank, height, Fraction(height, denominator))
    assert got == expected_row, ((m, d), got, expected_row)
    a_m = m * 2 ** (m - 2) - (m // 2) * comb(m, m // 2)
    assert height == Fraction(d * a_m, 2)

for a in range(2, 13):
    assert data([2, 2, 2 * a, 2 * a]) == (3, 2, 4)
    assert data([1, 1, a, a]) == (2, 2, 2)

print('Passed: all equal-degree dimension/rank/height/exponent table entries.')
print('Passed: both unequal four-row families for2<=a<=12.')
print('Rank proofs and arbitrary-degree family factorizations are in the companion note.')
