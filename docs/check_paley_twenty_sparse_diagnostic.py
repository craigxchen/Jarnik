"""Bounded sparse-certificate diagnostic for the Paley Hadamard matrix of order 20.

This checker covers supports 4 and 6 only. It is not a full-character
classification, nor evidence by itself for an unbounded-order theorem.
"""

from itertools import combinations


Q = 19
M = Q + 1


def paley_matrix():
    residues = {x * x % Q for x in range(1, Q)}

    def chi(x):
        x %= Q
        return 0 if x == 0 else (1 if x in residues else -1)

    return [[1] * M] + [
        [1] + [-(1 if i == j else chi(i - j)) for j in range(Q)]
        for i in range(Q)
    ]


def minimum_raw_transform(h, support_size):
    minimum = None
    count = 0
    for support in combinations(range(M), support_size):
        # Fix the first support entry positive to identify c and -c.
        for rest_positive in combinations(support[1:], support_size // 2 - 1):
            positive = set(rest_positive) | {support[0]}
            raw = [sum((1 if x in positive else -1) * h[x][a]
                       for x in support)
                   for a in range(1, M)]
            value = sum(map(abs, raw))
            minimum = value if minimum is None else min(minimum, value)
            count += 1
    return minimum, count


def main():
    h = paley_matrix()
    assert all(sum(h[x][a] * h[y][a] for a in range(M))
               == (M if x == y else 0)
               for x in range(M) for y in range(M))
    assert all(sum(h[x][a] for x in range(M)) == 0 for a in range(1, M))
    assert minimum_raw_transform(h, 4) == (24, 14_535)
    assert minimum_raw_transform(h, 6) == (28, 387_600)
    print("PASS: order-20 Paley; 14,535 support-4 and 387,600 support-6 "
          "sign-normalized certificates; raw transform minima 24 and 28.")


if __name__ == "__main__":
    main()
