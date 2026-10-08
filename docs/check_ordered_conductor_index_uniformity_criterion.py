"""Exact audits for the ordered conductor-index uniformity criterion."""

from fractions import Fraction
from itertools import combinations, product
from math import comb


def choose(n, k):
    return comb(n, k) if n >= k else 0


def f(t):
    return choose((t + 4) // 2, 5) + choose((t + 3) // 2, 5)


def subsequences(word, pattern):
    counts = [1] + [0] * len(pattern)
    for bit in word:
        for j in range(len(pattern) - 1, -1, -1):
            if bit == pattern[j]:
                counts[j + 1] += counts[j]
    return counts[-1]


def patterns(word):
    return subsequences(word, (0, 1, 1, 1, 0)) + subsequences(
        word, (1, 0, 0, 0, 1)
    )


def e4(signs):
    n, total = len(signs), sum(signs)
    numerator = total ** 4 - (6 * n - 8) * total ** 2 + 3 * n * (n - 2)
    assert numerator % 24 == 0
    return numerator // 24


def beta(m):
    return Fraction(1, 16) * (
        1 - Fraction(10, m - 1) - Fraction(5, m - 4)
        - Fraction(15, (m - 2) * (m - 3))
    )


def gamma(m):
    return Fraction(1, 16) * (
        1 - Fraction(10, m - 1)
        - Fraction(15, (m - 2) * (m - 3))
    )


def quartic_identity_checks():
    positive_pairs = ((1, 2), (1, 3), (2, 3), (0, 4))
    for word in product((0, 1), repeat=5):
        s = tuple(2 * x - 1 for x in word)
        quadratic = sum((1 if ij in positive_pairs else -1) * s[ij[0]] * s[ij[1]]
                        for ij in combinations(range(5), 2))
        quartic = sum((-1 if omitted in (0, 4) else 1)
                      * product_of(s[i] for i in range(5) if i != omitted)
                      for omitted in range(5))
        assert 16 * patterns(word) == 1 + quadratic + quartic
        all_pairs = sum(s[i] * s[j] for i, j in combinations(range(5), 2))
        assert 16 * int(len(set(word)) == 1) == 1 + all_pairs + e4(s)
    for n in range(5, 11):
        for word in product((0, 1), repeat=n):
            s = tuple(2 * x - 1 for x in word)
            direct = sum((2 * (indices[-1] - indices[0]) - n - 2)
                         * product_of(s[i] for i in indices)
                         for indices in combinations(range(n), 4))
            compressed = (n - 4) * e4(s) - 2 * sum(
                e4(s[:r]) + e4(s[r:]) for r in range(1, n)
            )
            assert direct == compressed


def product_of(values):
    result = 1
    for value in values:
        result *= value
    return result


def orthogonal_weighted_fixtures():
    for m in (5, 8, 16, 32, 34, 64):
        order = 1
        while order < m:
            order *= 2
        words = [tuple(bin(i & a).count('1') % 2 for i in range(m))
                 for a in range(1, order)]
        weight = len(words)
        assert all(sum((2 * word[i] - 1) * (2 * word[j] - 1) for word in words) == -1
                   for i, j in combinations(range(m), 2))
        actual = Fraction(sum(patterns(word) for word in words), comb(m, 5))
        assert actual >= beta(m) * weight
        constant = Fraction(sum(choose(sum(word), 5) + choose(m - sum(word), 5)
                                for word in words), comb(m, 5))
        assert constant >= gamma(m) * weight
        for delta in (Fraction(0), Fraction(1, 32), Fraction(1, 16), Fraction(1, 15)):
            assert actual + delta * constant >= (beta(m) + delta * gamma(m)) * weight


def main():
    checked = 0
    for n in range(1, 17):
        for word in product((0, 1), repeat=n):
            transitions = sum(x != y for x, y in zip(word, word[1:]))
            count = patterns(word)
            assert count >= f(transitions)
            if count == 0:
                assert transitions <= 5
            checked += 1
    for t in range(101):
        word = tuple(i % 2 for i in range(t + 1))
        assert patterns(word) == f(t)
        assert f(t + 1) - f(t) == choose((t + 3) // 2, 4)
        assert f(t + 2) - f(t + 1) >= f(t + 1) - f(t)
    quartic_identity_checks()
    orthogonal_weighted_fixtures()
    assert beta(33) < Fraction(1, 32) < beta(34)
    for m in range(34, 2001):
        assert beta(m) > Fraction(1, 32)
        assert beta(m + 1) > beta(m)
        assert gamma(m + 1) > gamma(m)
        assert Fraction(m * comb(m - 2, 3), 2 * comb(m, 5)) == Fraction(10, m - 1)
        assert Fraction(comb(m, 4), comb(m, 5)) == Fraction(5, m - 4)
        assert Fraction((m - 4) * m * (m - 1), 8 * comb(m, 5)) == Fraction(15, (m - 2) * (m - 3))
    for delta in (Fraction(0), Fraction(1, 32), Fraction(1, 16), Fraction(1, 16) + Fraction(1, 1000)):
        eta = (1 - 15 * delta) / 32
        assert eta > 0
        m = 34
        while beta(m) + delta * gamma(m) - delta < eta:
            m *= 2
        assert beta(m + 1) + delta * gamma(m + 1) - delta > eta
    print(f"PASS: {checked} binary words, alternating counts and convexity")
    print("PASS: quartic and constant-pattern identities, weighted fixtures, M=34 cutoff")
    print("PASS: intrinsic normalization improves the sufficient exponent to delta<1/15")
    print("PASS: ordered-sampling criterion; five-point arithmetic upper estimate remains unproved")


if __name__ == "__main__":
    main()
