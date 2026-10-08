"""Exact tableau checks for transpose symmetry and boundary comparisons.

The all-parameter proofs are in terminal_fusion_transpose_and_boundary_comparison.md.
This checker independently builds rectangular tableau paths rather than
calling the general Dynkin-coordinate fusion-grid implementation.
"""

from fractions import Fraction
from math import comb, factorial


def successors(shape, width):
    for row, length in enumerate(shape):
        if length < width and (row == 0 or shape[row - 1] > length):
            target = list(shape)
            target[row] += 1
            yield tuple(target)


def tableau_data(n, d, restricted):
    m = n * d
    start, end = (0,) * n, (d,) * n
    forward = [{start: 1}]
    for _ in range(m):
        nxt = {}
        for shape, count in forward[-1].items():
            for target in successors(shape, d):
                if restricted and target[0] - target[-1] >= d:
                    continue
                nxt[target] = nxt.get(target, 0) + count
        forward.append(nxt)
    rank = forward[m][end]
    backward = {end: 1}
    traces = [Fraction(0)] * (m + 1)
    for s in range(m, -1, -1):
        if s != m:
            backward = {
                shape: sum(backward.get(target, 0)
                           for target in successors(shape, d))
                for shape in forward[s]
            }
        assert sum(count * backward.get(shape, 0)
                   for shape, count in forward[s].items()) == rank
        traces[s] = sum(
            (Fraction(n * sum(length * (length + n - 1 - 2 * row)
                              for row, length in enumerate(shape)) - s * s, n)
             * count * backward.get(shape, 0)
             for shape, count in forward[s].items()), Fraction(0))
    return rank, traces


def source_dimension(n, d):
    numerator = factorial(n * d)
    denominator = 1
    for i in range(n):
        numerator *= factorial(i)
        denominator *= factorial(d + i)
    assert numerator % denominator == 0
    return numerator // denominator


def case(n, d):
    m, kappa = n * d, n + d - 1
    full_rank, full_traces = tableau_data(n, d, False)
    h, traces = tableau_data(n, d, True)
    assert full_rank == source_dimension(n, d)
    fundamental = Fraction(n * n - 1, n)
    means = [Fraction(s * (m - s), m - 1) * fundamental
             for s in range(m + 1)]
    assert full_traces == [full_rank * value for value in means]
    deficits = [h * means[s] - traces[s] for s in range(m + 1)]
    assert deficits == list(reversed(deficits))
    delta = (m - 1) * deficits[2] / (2 * kappa)
    orders = [(comb(s, 2) * deficits[2] - deficits[s]) / (2 * kappa)
              for s in range(m + 1)]
    assert delta.denominator == 1 and delta > 0
    assert all(value.denominator == 1 and value >= 0 for value in orders)
    B = m * delta * 2 ** (m - 3) - sum(
        comb(m, s) * orders[s] for s in range(m + 1))
    assert B == sum(comb(m, s) * deficits[s]
                    for s in range(m + 1)) / (2 * kappa)
    return dict(h=h, full_rank=full_rank, J=deficits,
                delta=delta, orders=orders, B=B)


def boundary_numerator(j):
    return Fraction(j * 4 ** (j - 1)) - Fraction(j * comb(2 * j, j), 2)


def threshold(n, d):
    return 2 ** (n * d) - 2 * (2 ** d - 1) ** n + (2 ** d - 2) ** n


def main():
    data = {(n, d): case(n, d)
            for n in range(2, 7) for d in range(2, 7)}
    for (n, d), first in data.items():
        second = data[d, n]
        assert first['h'] + second['h'] == first['full_rank']
        for key in ('J', 'delta', 'orders', 'B'):
            assert first[key] == second[key]
        if d == 2:
            assert first['h'] == first['delta'] == 1
            assert first['orders'] == [max(0, s - n)
                                       for s in range(2 * n + 1)]
            assert first['B'] == boundary_numerator(n)
        if n == 2:
            assert first['h'] == comb(2 * d, d) // (d + 1) - 1
            assert first['B'] == boundary_numerator(d)
    assert [boundary_numerator(n) - threshold(n, 2)
            for n in range(2, 8)] == [0, 0, 6, 80, 670, 4522]
    equalities = set()
    for j in range(2, 101):
        central = comb(2 * j, j)
        assert central * central * (3 * j + 1) <= 16 ** j
        assert (4 * (j + 1) ** 2 * (3 * j + 1)
                - (2 * j + 1) ** 2 * (3 * j + 4)) == j
        B = boundary_numerator(j)
        h_two_rows = central // (j + 1) - 1
        assert B >= 2 * h_two_rows
        assert B >= threshold(j, 2)
        if B == 2 * h_two_rows:
            equalities.add((2, j))
        if B == threshold(j, 2):
            equalities.add((j, 2))
    assert equalities == {(2, 2), (3, 2)}
    print('PASS: 25 tableau grids verify all Casimir-deficit and transpose identities.')
    print('PASS: boundary formulas through size 100; equalities only (2,2),(3,2).')


if __name__ == '__main__':
    main()
