#!/usr/bin/env python3
"""Exact graph check for the all-degree homogeneous-binomial bound."""

from itertools import combinations, product


EDGES = list(combinations(range(8), 2))
CUTS = list(range(1, 128))  # one representative of each cut modulo complement
INSIDE = [[k for k, (i, j) in enumerate(EDGES) if s >> i & 1 and s >> j & 1]
          for s in CUTS]
CROSS = [sum(1 << k for k, (i, j) in enumerate(EDGES) if ((s >> i) ^ (s >> j)) & 1)
         for s in CUTS]


def canonical(cut):
    return min(cut, cut ^ 255)


def check_matrix(block):
    matrix = [[0] * 4 for _ in range(4)]
    for i in range(3):
        for j in range(3):
            matrix[i][j] = block[3 * i + j]
    for i in range(3):
        matrix[i][3] = -sum(matrix[i][:3])
    for j in range(3):
        matrix[3][j] = -sum(matrix[i][j] for i in range(3))
    matrix[3][3] = sum(block)
    assert all(sum(row) == 0 for row in matrix)
    assert all(sum(matrix[i][j] for i in range(4)) == 0 for j in range(4))
    edges = [matrix[i][j - 4] if i < 4 <= j else 0 for i, j in EDGES]
    support = sum(1 << k for k, value in enumerate(edges) if value)
    assert support
    orders = [sum(edges[k] for k in inside) for inside in INSIDE]
    for s, value in zip(CUTS, orders):
        complement_value = sum(b for (i, j), b in zip(EDGES, edges)
                               if not (s >> i & 1) and not (s >> j & 1))
        assert value == complement_value
    good = [s for s, crossing in zip(CUTS, CROSS) if support & ~crossing == 0]
    assert good  # the defining four-four bipartition is good
    neighbors = {i: [] for i in range(8)}
    for (i, j), value in zip(EDGES, edges):
        if value:
            neighbors[i].append(j)
            neighbors[j].append(i)
    assert all(len(neighbors[i]) != 1 for i in range(8))
    j = next(i for i in range(8) if neighbors[i])
    i, k = neighbors[j][:2]
    first = {canonical(s ^ (1 << i) ^ (1 << j)) for s in good}
    second = {canonical(s ^ (1 << j) ^ (1 << k)) for s in good}
    assert len(first) == len(second) == len(good)
    assert not first.intersection(second)
    b_ij = edges[EDGES.index(tuple(sorted((i, j))))]
    b_jk = edges[EDGES.index(tuple(sorted((j, k))))]
    assert all(orders[s - 1] == -b_ij for s in first)
    assert all(orders[s - 1] == -b_jk for s in second)
    absolute_sum = sum(map(abs, orders))
    assert absolute_sum >= (abs(b_ij) + abs(b_jk)) * len(good) >= 2 * len(good)
    return len(good), absolute_sum == 2 * len(good)


if __name__ == "__main__":
    checked = equalities = 0
    good_counts = set()
    for block in product((-1, 0, 1), repeat=9):
        if not any(block):
            continue
        count, equality = check_matrix(block)
        checked += 1
        equalities += equality
        good_counts.add(count)
    assert checked == 19682
    print(f"PASS: {checked} exact signed matrices; complementary orders and both disjoint injections.")
    print(f"PASS: cancellation-loss inequality; {equalities} equality cases; good-cut counts {sorted(good_counts)}.")
