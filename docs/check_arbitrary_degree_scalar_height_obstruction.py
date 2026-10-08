"""Exact finite certificates for the arbitrary-degree scalar obstruction.

The companion proof derives the cone rays and the complete 32-vector
cutoff. This script does not replace the general reduction proof.
"""

from collections import Counter
from fractions import Fraction
from math import comb


def subset_counts(degrees):
    counts = Counter({Fraction(0): 1})
    for degree in degrees:
        counts += Counter({weight + degree: n for weight, n in counts.items()})
    assert sum(counts.values()) == 2 ** len(degrees)
    return counts


def height(degrees):
    degrees = tuple(map(Fraction, degrees))
    total = sum(degrees)
    counts = subset_counts(degrees)
    value = total * 2 ** (len(degrees) - 3) - sum(
        n * max(0, weight - total / 2) for weight, n in counts.items()
    )
    absolute_sum = sum(n * abs(2 * weight - total)
                       for weight, n in counts.items())
    assert 4 * value == total * 2 ** (len(degrees) - 1) - absolute_sum
    return value


def dimension(degrees):
    total = sum(degrees)
    assert total % 2 == 0
    coefficients = [1]
    for degree in degrees:
        # Sliding window computes multiplication by 1+x+...+x^degree.
        nxt, window = [], 0
        for i in range(len(coefficients) + degree):
            if i < len(coefficients):
                window += coefficients[i]
            if i - degree - 1 >= 0:
                window -= coefficients[i - degree - 1]
            nxt.append(window)
        coefficients = nxt
    middle = total // 2
    return coefficients[middle] - coefficients[middle - 1]


def data(degrees):
    total = sum(degrees)
    assert total % 2 == 0
    counts = subset_counts(degrees)
    assert counts[total // 2] % 2 == 0
    return dimension(degrees), counts[total // 2] // 2, height(degrees)


def rays(m):
    result = []
    for q in range(4, m + 1):
        prefix = [Fraction(0)] * (m - q)
        result.append(prefix + [Fraction(1)] * q)
        if q > 4:
            result.extend([
                prefix + [Fraction(1)] * (q - 3) + [Fraction(q - 3)] * 3,
                prefix + [Fraction(1)] * (q - 2) + [Fraction(q - 2, 2)] * 2,
                prefix + [Fraction(1)] * (q - 1) + [Fraction(q - 3)],
            ])
    for vector in result:
        assert vector == sorted(vector)
        assert 2 * (vector[-1] + vector[-2]) <= sum(vector)
    return result


def partitions(total, length, least=1):
    if length == 1:
        if total >= least:
            yield (total,)
        return
    for first in range(least, total // length + 1):
        for remaining in partitions(total - first, length - 1, first):
            yield (first,) + remaining


def reduce_degrees(degrees):
    degrees = tuple(sorted(degrees))
    while degrees[-1] + degrees[-2] > sum(degrees) // 2:
        middle = sum(degrees) // 2
        assert degrees[-1] < middle
        defect = degrees[-1] + degrees[-2] - middle
        before = data(degrees)
        reduced = list(degrees)
        reduced[-1] -= defect
        reduced[-2] -= defect
        degrees = tuple(sorted(reduced))
        assert min(degrees) > 0
        assert max(degrees) < sum(degrees) // 2
        after = data(degrees)
        assert after == (before[0], before[1] + 1, before[2])
    return degrees


EXPECTED = {
    5: (Fraction(1), Fraction(10), 4, Fraction(8, 3)),
    6: (Fraction(7, 4), Fraction(80, 7), 5, Fraction(4)),
    7: (Fraction(3), Fraction(35, 3), 4, Fraction(29, 6)),
    8: (Fraction(31, 6), Fraction(420, 31), 8, Fraction(132, 19)),
    9: (Fraction(9), Fraction(14), 11, Fraction(8)),
}


def main():
    total_checked = 0
    for m, (expected_ratio, expected_bound, expected_count, expected_minimum) in EXPECTED.items():
        ray_ratio = min(height(vector) / sum(vector) for vector in rays(m))
        assert ray_ratio == expected_ratio
        degree_bound = Fraction(comb(m, m // 2), 1) / ray_ratio
        assert degree_bound == expected_bound
        count, minima = 0, []
        maximum_total = degree_bound.numerator // degree_bound.denominator
        for total in range(m + m % 2, maximum_total + 1, 2):
            for vector in partitions(total, m):
                if 2 * (vector[-1] + vector[-2]) > total:
                    continue
                rank, cuts, value = data(vector)
                denominator = min(cuts, rank - 1)
                assert value > 2 * denominator, (vector, rank, cuts, value)
                assert value >= ray_ratio * total
                if denominator:
                    minima.append(value / denominator)
                else:
                    assert vector == (2, 2, 2, 2, 2)
                    assert (rank, cuts, value) == (6, 0, 10)
                count += 1
        assert count == expected_count, (m, count)
        assert min(minima) == expected_minimum
        total_checked += count
    assert total_checked == 32

    # Finite ray formulas appearing in the large-row proof.
    for q in range(4, 13):
        small_abs = lambda n: Fraction(sum(comb(n, j) * abs(n - 2 * j)
                                         for j in range(n + 1)), 2**n)
        full = [vector for vector in rays(q) if vector[0] == 1]
        expected = [2 ** (q - 3) * (q - 2 * small_abs(q))]
        if q > 4:
            expected += [(q - 3) * 2 ** (q - 3),
                         2 ** (q - 3) * (q - 2 - small_abs(q - 2)),
                         2 ** (q - 2) - 2]
        assert [height(vector) for vector in full] == expected
        if q >= 6:
            assert min(expected) == 2 ** (q - 2) - 2

    # Nontrivial inverse bracket removals with arbitrarily chosen large
    # coefficients; these are checks of the formulas, not search bounds.
    forced_examples = 0
    for m in range(4, 10):
        canonical = [1] * (m - 1) + [m - 3]
        for increment in (1, 7, 1000):
            lifted = canonical.copy()
            lifted[-2] += increment
            lifted[-1] += increment
            assert reduce_degrees(lifted) == tuple(canonical)
            forced_examples += 1
    for vector in ((1, 5, 8, 12), (2, 4, 6, 8), (2, 3, 7, 8)):
        residual = reduce_degrees(vector)
        assert len(set(residual)) == 1
        rank, cuts, value = data(vector)
        assert value == 2 * (rank - 1)
        assert value >= 2 * min(cuts, rank - 1)
        forced_examples += 1
    for c in range(1, 9):
        assert data((c, c, c, c)) == (c + 1, 3, 2 * c)

    print("Passed: all 32 residual degree vectors and all five cone ratio minima.")
    print("Passed: finite ray formulas, four-row identities, and "
          f"{forced_examples} forced-factor examples.")
    print("The cone decomposition and the m >= 10 proof are in the companion note.")


if __name__ == "__main__":
    main()
