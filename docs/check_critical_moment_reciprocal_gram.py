"""Exact reciprocal-Gram identities; no full eight-row solution is asserted."""

from fractions import Fraction as F
from itertools import combinations


def determinant(matrix):
    a = [list(map(F, row)) for row in matrix]
    result = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result = -result
        scale = a[i][i]
        result *= scale
        for j in range(i + 1, len(a)):
            factor = a[j][i] / scale
            for k in range(i + 1, len(a)):
                a[j][k] -= factor * a[i][k]
    return result


def vandermonde(values):
    result = F(1)
    for i, j in combinations(range(len(values)), 2):
        result *= values[j] - values[i]
    return result


def polynomial(roots):
    out = [F(1)]
    for root in roots:
        new = [F(0)] * (len(out) + 1)
        for i, value in enumerate(out):
            new[i] -= root * value
            new[i + 1] += value
        out = new
    return out


def lorentz_distance(x, y):
    delta = [F(a) - F(b) for a, b in zip(x, y)]
    return sum(t * t for t in delta[:-1]) - delta[-1] ** 2


def run():
    a = list(map(F, [0, 1, 2, 3]))
    b = list(map(F, [-2, -1, 4, 5]))
    gamma = F(7)
    c = [[gamma - (x - y) ** 2 / 2 for y in b] for x in a]
    k = [[gamma, 0, F(-1, 2)], [0, 1, 0], [F(-1, 2), 0, 0]]
    assert determinant(k) == F(-1, 4)
    assert determinant(c) == 0
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            va, vb = [1, x, x * x], [1, y, y * y]
            assert c[i][j] == sum(va[r] * k[r][s] * vb[s]
                                  for r in range(3) for s in range(3))
    for rows in combinations(range(4), 3):
        for cols in combinations(range(4), 3):
            minor = determinant([[c[i][j] for j in cols] for i in rows])
            assert minor == -vandermonde([a[i] for i in rows]) * vandermonde(
                [b[j] for j in cols]) / 4
            assert minor < 0
    row_means = [sum(row) / 4 for row in c]
    col_means = [sum(c[i][j] for i in range(4)) / 4 for j in range(4)]
    mean = sum(row_means) / 4
    for i in range(4):
        for j in range(4):
            assert c[i][j] - row_means[i] - col_means[j] + mean == (
                a[i] - sum(a) / 4) * (b[j] - sum(b) / 4)

    eta_a = list(map(F, [2, -3, 7, 1]))
    eta_b = list(map(F, [-1, 4, 2, 8]))
    gamma4 = F(11)
    fourth = [[gamma4 - F(2, 3) * (x - y) * (eta_a[i] - eta_b[j])
               + (x - y) ** 4 / 24 for j, y in enumerate(b)]
              for i, x in enumerate(a)]
    bordered = [row + [F(1), a[i]] for i, row in enumerate(fourth)] + [
        [F(1)] * 4 + [F(0)] * 2, b + [F(0)] * 2]
    assert determinant(bordered) == 0
    for rows in combinations(range(4), 3):
        wa = {i: 1 / product(a[i] - a[j] for j in rows if j != i) for i in rows}
        assert sum(wa.values()) == 0
        assert sum(wa[i] * a[i] for i in rows) == 0
        assert sum(wa[i] * a[i] ** 2 for i in rows) == 1
        for cols in combinations(range(4), 3):
            wb = {i: 1 / product(b[i] - b[j] for j in cols if j != i) for i in cols}
            assert sum(wa[i] * fourth[i][j] * wb[j]
                       for i in rows for j in cols) == F(1, 4)
            row_ids, col_ids = list(rows) + [4, 5], list(cols) + [4, 5]
            minor = determinant([[bordered[i][j] for j in col_ids] for i in row_ids])
            assert minor == vandermonde([a[i] for i in rows]) * vandermonde(
                [b[j] for j in cols]) / 4
            assert minor > 0

    for roots, kind in [([-3, 1, 2], 0), ([-4, -1, 2, 3], 1),
                        ([-6, 1, 2, 3], -1)]:
        coeff = polynomial(roots)
        p1 = sum(F(1, z) for z in roots)
        p2 = sum(F(1, z * z) for z in roots)
        p3 = sum(F(1, z ** 3) for z in roots)
        assert coeff[1] / coeff[0] == -p1
        assert 2 * coeff[2] / coeff[0] == p1 * p1 - p2
        if kind == 0:
            assert coeff[2] == 0 and p1 * p1 == p2
        else:
            assert all(coeff[j] == 0 for j in range(3, len(coeff), 2))
            assert p1 ** 3 - 3 * p1 * p2 + 2 * p3 == 0
            delta_beta, delta_eta = 2 * p1, 2 * p3
            assert delta_eta == -delta_beta ** 3 / 8 + F(3, 4) * delta_beta * (2 * p2)
            distance = 4 * (p2 - p1 * p1)
            assert distance == -8 * coeff[2] / coeff[0]
            assert kind * distance > 0

    for roots in [[-3, 1, 2], [1, 5, 9, -7, -8]]:
        p1, p2, p3, p4 = [sum(F(1, z ** k) for z in roots) for k in range(1, 5)]
        assert p1 ** 2 == p2
        assert 4 * p1 * p3 == p2 ** 2 + 3 * p4
        delta_beta, delta_eta = 2 * p1, 2 * p3
        # No common roots: gamma4=p4 and C4=-p4 for these pair fixtures.
        assert -p4 == p4 - F(2, 3) * delta_beta * delta_eta + delta_beta ** 4 / 24

    c_order = [6, 7, 0, 1, 2, 3, 4, 5]
    beta_order = [4, 5, 0, 1, 2, 3, 6, 7]
    c_rank = {i: r for r, i in enumerate(c_order)}
    beta_rank = {i: r for r, i in enumerate(beta_order)}
    for i, j in combinations(range(8), 2):
        cross = (i < 4) != (j < 4)
        successor_b = i >= 4 and (i < 6) != (j < 6)
        multiplier = -1 if cross or successor_b else 1
        assert (beta_rank[i] - beta_rank[j]) * multiplier * (
            c_rank[i] - c_rank[j]) > 0

    left = [(1, 0, 0, 0, 0), (-1, 0, 0, 0, 0),
            (0, 1, 0, 0, 0), (0, 0, 1, 0, 0)]
    right = [(0, 0, 0, 0, 1), (0, 0, 0, F(3, 4), F(5, 4)),
             (0, 0, 0, 0, -1), (0, 0, 0, F(3, 4), F(-5, 4))]
    assert all(lorentz_distance(x, y) == 0 for x in left for y in right)
    assert all(lorentz_distance(left[i], left[j]) > 0
               for i, j in combinations(range(4), 2))
    for i, j in combinations(range(4), 2):
        sign = 1 if (i < 2) == (j < 2) else -1
        assert sign * lorentz_distance(right[i], right[j]) > 0
    print("PASS: exact reciprocal Gram factorization, 16 minors, and centering")
    print("PASS: fourth-order divided differences, 16 bordered minors, and within-group identity")
    print("PASS: three pair residuals, reciprocal order, and causal-sign fixture")


def product(values):
    out = F(1)
    for value in values:
        out *= value
    return out


if __name__ == "__main__":
    run()
