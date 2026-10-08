"""Bounded audits for the Paley square-root-support fourth-moment proof.

The asymptotic support theorem uses Weil's character-sum bound; these
fixtures check its moment normalization and the norm-one trade dictionary.
"""

from itertools import combinations
from math import prod, sqrt


def paley_matrix(q):
    residues = {x * x % q for x in range(1, q)}

    def chi(z):
        z %= q
        return 0 if z == 0 else (1 if z in residues else -1)

    return [[1] * (q + 1)] + [
        [1] + [-(1 if x == a else chi(x - a)) for a in range(q)]
        for x in range(q)
    ]


def hadamard(h):
    m = len(h)
    return all(sum(h[x][a] * h[y][a] for a in range(m))
               == (m if x == y else 0)
               for x in range(m) for y in range(m))


def transform(h, c):
    m = len(h)
    return [sum(c[x] * h[x][a] for x in range(m)) for a in range(m)]


def moment_fixture(h, support_size):
    m = len(h)
    support = tuple(range(support_size))
    c = [0] * m
    for x in support[:support_size // 2]:
        c[x] = 1
    for x in support[support_size // 2:]:
        c[x] = -1
    t = transform(h, c)
    assert t[0] == 0
    assert sum(z * z for z in t) == m * support_size

    quartic = sum(
        prod(c[x] for x in rows)
        * sum(prod(h[x][a] for x in rows) for a in range(m))
        for rows in combinations(support, 4)
    )
    fourth = sum(z ** 4 for z in t)
    assert fourth == m * (3 * support_size ** 2 - 2 * support_size) + 24 * quartic

    # The L1--L2--L4 interpolation is exact after clearing denominators.
    absolute = sum(abs(z) for z in t)
    assert support_size ** 3 * m ** 3 <= absolute ** 2 * fourth
    return support_size, absolute, fourth


def norm_one_trade_fixture():
    q = 11
    h = paley_matrix(q)
    m = q + 1
    a, b = 1, 2  # Two nonconstant column labels.
    c = [(h[x][a] + h[x][b]) // 2 for x in range(m)]
    assert sum(c) == 0 and set(c) <= {-1, 0, 1}
    support = [x for x in range(m) if c[x]]
    t = transform(h, c)
    active = [j for j in range(m) if t[j]]
    assert len(support) == 6 and active == [a, b]
    assert all(abs(t[j]) == len(support) for j in active)
    assert len(support) * len(active) == m
    assert all(h[x][j] == c[x] * (1 if t[j] > 0 else -1)
               for x in support for j in active)

    # Switching this rank-one rectangle gives another Hadamard matrix.
    switched = [row[:] for row in h]
    for x in support:
        for j in active:
            switched[x][j] *= -1
    assert hadamard(switched)

    copies = 5
    physical = [[h[x][j] for j in range(1, m) for _ in range(copies)]
                for x in range(m)]
    used = [0] * m
    assigned = []
    for x in range(m):
        label = 1 + x % q
        j = copies * (label - 1) + used[label]
        used[label] += 1
        physical[x][j] *= -1
        assigned.append(label)
    assert max(used) == 2
    coefficients = [sum(c[x] * physical[x][j] for x in range(m)) // 2
                    for j in range(copies * q)]
    assert all(sum(c[x] * physical[x][j] for x in range(m)) % 2 == 0
               for j in range(copies * q))
    matching = sum(x in support and assigned[x] in active for x in range(m))
    coefficient_sum = sum(abs(z) for z in coefficients)
    assert coefficient_sum == copies * m // 2 + len(support) - 2 * matching
    assert coefficient_sum > copies * q / 2
    return len(support), len(active), matching, coefficient_sum


def explicit_cutoff_fixture():
    q = 4099
    h = 4
    assert q % 4 == 3 and all(q % d for d in range(2, int(sqrt(q)) + 1))
    assert h <= sqrt(q) / 16
    delta = (3 * sqrt(q) + 5) / (q + 1)
    assert delta <= 4 / sqrt(q)
    assert (3 / h - 2 / h ** 2 + delta * h) < 7 / 8
    assert sqrt(8 / 7) > 17 / 16
    return q


def main():
    fixtures = 0
    for q in (11, 19, 31, 43):
        h = paley_matrix(q)
        assert hadamard(h)
        for support_size in (4, 6, 8):
            moment_fixture(h, support_size)
            fixtures += 1
    trade = norm_one_trade_fixture()
    cutoff = explicit_cutoff_fixture()
    print("PASS: %d Paley moment fixtures; order-12 norm-one %d-by-%d "
          "rectangle, matching count %d, B=%d; explicit q=%d cutoff."
          % (fixtures, *trade, cutoff))


if __name__ == "__main__":
    main()
