"""Exhaustive integer certificate audit for the order-twelve Paley core."""

from collections import Counter
from itertools import combinations


Q = 11
M = 12
B = 5


def paley_matrix():
    residues = {x * x % Q for x in range(1, Q)}

    def chi(x):
        x %= Q
        return 0 if not x else (1 if x in residues else -1)

    return [[1] * M] + [
        [1] + [-(1 if i == j else chi(i - j)) for j in range(Q)]
        for i in range(Q)
    ]


def normalize(c):
    sign = next(x for x in c if x)
    assert sign in (-1, 1)
    return tuple(sign * x for x in c)


def expected_norm_one(h):
    expected = set()
    for x, y in combinations(range(M), 2):
        c = [0] * M
        c[x], c[y] = 1, -1
        expected.add(normalize(c))
    columns = [tuple(h[x][a] for x in range(M)) for a in range(1, M)]
    expected.update(normalize(c) for c in columns)
    for a, b in combinations(columns, 2):
        for sign in (-1, 1):
            expected.add(normalize(tuple((x + sign * y) // 2
                                         for x, y in zip(a, b))))
    assert len(expected) == 187
    return expected


def physical_height(c, h, labels):
    assigned = {a: [x for x in range(M) if labels[x] == a]
                for a in range(1, M)}
    assert max(map(len, assigned.values())) <= 2
    total = 0
    for a in range(1, M):
        twice = sum(c[x] * h[x][a] for x in range(M))
        assert twice % 2 == 0
        u = twice // 2
        total += (B - len(assigned[a])) * abs(u)
        total += sum(abs(u - c[x] * h[x][a]) for x in assigned[a])
    return total


def main():
    h = paley_matrix()
    assert all(sum(h[x][a] * h[y][a] for a in range(M))
               == (M if x == y else 0)
               for x in range(M) for y in range(M))
    assert all(sum(h[x][a] for x in range(M)) == 0 for a in range(1, M))
    actual_one = set()
    table = {}
    seen = 0
    for pairs in range(1, M // 2 + 1):
        norms = Counter()
        for support in combinations(range(M), 2 * pairs):
            for rest_positive in combinations(support[1:], pairs - 1):
                positive = set(rest_positive) | {support[0]}
                c = tuple((1 if x in positive else -1) if x in support else 0
                          for x in range(M))
                assert sum(c) == 0 and normalize(c) == c
                raw = [sum(c[x] * h[x][a] for x in support) for a in range(M)]
                norm = sum(map(abs, raw))
                norms[norm] += 1
                seen += 1
                if norm == M:
                    actual_one.add(c)
                else:
                    assert B * norm - 2 * len(support) >= 2 * 36
        table[2 * pairs] = (min(norms), norms[M],
                            min([v for v in norms if v > M] or [0]))
    assert seen == 36894
    assert table == {2: (12, 66, 0), 4: (16, 0, 16), 6: (12, 110, 20),
                     8: (20, 0, 20), 10: (20, 0, 20), 12: (12, 11, 28)}
    assert actual_one == expected_norm_one(h)

    # Every listed equality case has aligned signs at all selected rows
    # for each of its nonzero transform labels; verify the proof of (3).
    for c in actual_one:
        raw = [sum(c[x] * h[x][a] for x in range(M)) for a in range(M)]
        active = [a for a in range(1, M) if raw[a]]
        support = [x for x in range(M) if c[x]]
        for x in support:
            for a in active:
                assert c[x] * h[x][a] * raw[a] > 0
        k_bound = min(len(support), 2 * len(active))
        assert 30 + len(support) - 2 * k_bound >= 28

    labels = tuple(1 + x % 11 for x in range(M))
    pair = tuple(1 if x == 1 else (-1 if x == 2 else 0) for x in range(M))
    assert physical_height(pair, h, labels) == 28
    assert 2 * 28 > B * (M - 1)
    print("PASS: 36,894 sign-normalized ternary certificates; exact support "
          "table; all 187 norm-one vectors and their flip signs; "
          "assignment-independent lower bound 28 > 55/2, attained canonically.")


if __name__ == "__main__":
    main()
