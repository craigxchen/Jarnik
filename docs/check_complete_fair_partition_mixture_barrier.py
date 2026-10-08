"""Exact balanced-cut means and the universal pair-mixture flip decrement."""

from fractions import Fraction
from itertools import combinations


def main():
    cases = 0
    for m in range(4, 17, 2):
        k = m // 2
        cuts = list(combinations(range(m), k))
        means = []
        for pairs in range(1, k + 1):
            total = 0
            for cut in cuts:
                plus = sum(x < pairs for x in cut)
                minus = sum(pairs <= x < 2 * pairs for x in cut)
                total += abs(plus - minus)
            means.append(Fraction(total, len(cuts)))
            cases += 1
        a2 = Fraction(k, 2 * k - 1)
        a4 = Fraction(k * (3 * k - 5), (2 * k - 1) * (2 * k - 3))
        assert means[0] == a2
        assert means[1] == a4
        assert all(x <= y for x, y in zip(means, means[1:]))
        assert a4 / a2 - 1 == Fraction(k - 2, 2 * k - 3)
        assert Fraction(1, 2) / a2 == Fraction(m - 1, m)

        pairs = list(combinations(range(m), 2))
        for selected in (set(range(k - 1)), set(range(k + 1))):
            separation = Fraction(sum((x in selected) != (y in selected)
                                      for x, y in pairs), len(pairs))
            decrement = Fraction(m - 1, m) * (a2 - separation)
            assert decrement == Fraction(2, m * m)
    print("PASS: %d exact balanced-certificate means for M=4,...,16; "
          "monotonicity, mixture constants, and both flip directions." % cases)


if __name__ == "__main__":
    main()
