"""Exact Euler-word generic-rank and Boolean specialization checks."""

from fractions import Fraction
from itertools import combinations


def determinant(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        lead = a[j][j]
        result *= lead
        for i in range(j + 1, len(a)):
            coeff = a[i][j] / lead
            a[i] = [x - coeff * y for x, y in zip(a[i], a[j])]
    return result


def rank(matrix):
    rows = [list(map(Fraction, row)) for row in matrix]
    pivots = 0
    for j in range(len(rows[0])):
        q = next((i for i in range(pivots, len(rows)) if rows[i][j]), None)
        if q is None:
            continue
        rows[pivots], rows[q] = rows[q], rows[pivots]
        lead = rows[pivots][j]
        rows[pivots] = [x / lead for x in rows[pivots]]
        for i in range(len(rows)):
            if i != pivots:
                coeff = rows[i][j]
                rows[i] = [x - coeff * y for x, y in zip(rows[i], rows[pivots])]
        pivots += 1
    return pivots


def bit_matrix(columns):
    return [[(v >> i) & 1 for v in columns] for i in range(7)]


def check_partition():
    orbit, v = [], 1
    while v not in orbit:
        orbit.append(v)
        v = ((v << 1) & 127) ^ (3 if v & 64 else 0)
    assert v == 1 and sorted(orbit) == list(range(1, 128))
    assert orbit[:7] == [1, 2, 4, 8, 16, 32, 64]
    for q in range(127):
        cols = [orbit[(q + j) % 127] for j in range(7)]
        det = determinant(bit_matrix(cols))
        assert det.denominator == 1 and int(det) % 2
    groups = [orbit[j:j + 7] for j in range(0, 127, 7)]
    assert len(groups) == 19
    assert all(rank(bit_matrix(group)) == len(group) for group in groups)
    slope = 20_805_120 - 19 * 806_400
    constant = 2 - 19 * 10
    assert (slope, constant) == (5_483_520, -188)
    assert slope > 0 and slope + constant > 0


def check_boolean_fixture():
    C = [[1, 0, 1, 0, 1, 0], [0, 1, 0, 0, 1, 1], [0, 0, 1, 1, 0, 1]]
    a, lam = [1, 3, 4, 5, 6, 9], [-32, 15, 56, -65, -24, 9]
    x, cols = [v * v for v in a], list(zip(*C))
    assert len(set(cols)) == 6 and all(any(col) for col in cols)
    for indices, expected in [([0, 1, 3], 1), ([2, 4, 5], 2)]:
        assert determinant([[row[j] for j in indices] for row in C]) == expected

    def jacobian(nodes):
        return C + [[b * v for b, v in zip(row, nodes)] for row in C]

    assert all(sum(b * v for b, v in zip(row, lam)) == 0 for row in jacobian(x))
    assert rank(jacobian(x)) == 5
    assert rank(jacobian([1, 2, 3, 4, 5, 6])) == 6
    assert determinant(jacobian([1, 2, 3, 4, 5, 6])) == 4
    assert len({Fraction(abs(v), t) for v, t in zip(lam, a)}) > 1
    count = 0
    for y in combinations(range(1, 10), 6):
        formula = -((y[2] - y[0]) * (y[4] - y[1]) * (y[5] - y[3])
                    + (y[4] - y[0]) * (y[5] - y[1]) * (y[2] - y[3]))
        assert determinant(jacobian(y)) == formula
        count += 1
    assert count == 84
    return count


if __name__ == "__main__":
    check_partition()
    count = check_boolean_fixture()
    print("Passed 127 binary bases, 19-group partition, and all-K degree budget.")
    print(f"Passed Boolean square-node kernel, exact ranks, and {count} determinant audits.")
