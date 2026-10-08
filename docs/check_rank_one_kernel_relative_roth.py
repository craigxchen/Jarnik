"""Exact finite checks for rank_one_kernel_relative_roth_obstruction.md.

Run from the repository root with

    python3 docs/check_rank_one_kernel_relative_roth.py

The proof's Diophantine separation is symbolic prose.  This checker verifies
the finite Paley-core linear algebra, all pair isolators, the derivative
constants of the four quadratic forms, and the cross-class matchings for all
integer partitions through 24 coordinates.
"""

from fractions import Fraction


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def paley_core_11():
    squares = {1, 3, 4, 5, 9}

    def chi(value):
        value %= 11
        if value == 0:
            return 0
        return 1 if value in squares else -1

    return [
        [-1 if x == y else chi(x - y) for y in range(11)]
        for x in range(11)
    ]


def check_paley_core():
    s = paley_core_11()
    st = transpose(s)
    identity = [[int(i == j) for j in range(11)] for i in range(11)]
    jmat = [[1 for _ in range(11)] for _ in range(11)]
    gram = [
        [12 * identity[i][j] - jmat[i][j] for j in range(11)]
        for i in range(11)
    ]

    assert [sum(row) for row in s] == [-1] * 11
    assert [sum(row) for row in st] == [-1] * 11
    assert matmul(s, st) == gram
    assert matmul(st, s) == gram

    inverse_numerator = [
        [st[i][j] - 1 for j in range(11)] for i in range(11)
    ]
    scaled_inverse_check = matmul(s, inverse_numerator)
    assert scaled_inverse_check == [
        [12 * identity[i][j] for j in range(11)] for i in range(11)
    ]

    assert all(row.count(1) == 5 and row.count(-1) == 6 for row in s)
    for a in range(11):
        for b in range(a):
            column_difference = [s[i][a] - s[i][b] for i in range(11)]
            assert set(column_difference) <= {-2, 0, 2}
            assert sum(value != 0 for value in column_difference) == 6
            assert sum(abs(Fraction(value, 12)) for value in column_difference) == 1

            isolated = [
                sum(Fraction(column_difference[i], 12) * s[i][j] for i in range(11))
                for j in range(11)
            ]
            expected = [Fraction(int(j == a) - int(j == b)) for j in range(11)]
            assert isolated == expected
            assert sum(Fraction(value, 12) for value in column_difference) == 0

    for a in range(11):
        for b in range(a):
            assert sum(s[a][j] != s[b][j] for j in range(11)) == 6

    return 55


def check_quadratic_forms():
    # Q(x,y)=a*x^2+b*x*y+c*y^2 restricts on the unit circle with
    # derivative amplitude sqrt(b^2+(c-a)^2).
    forms = {
        "negative_steep_or_shallow": (1, 4, 1),
        "positive_steep_or_shallow": (1, -4, 1),
        "slope_one_over_sqrt_three": (1, 0, -3),
        "slope_sqrt_three": (3, 0, -1),
    }
    for a, b, c in forms.values():
        assert b * b + (c - a) * (c - a) <= 16

    irrational_target_forms = {
        1: forms["positive_steep_or_shallow"],
        2: forms["slope_one_over_sqrt_three"],
        4: forms["slope_sqrt_three"],
        5: forms["positive_steep_or_shallow"],
        7: forms["negative_steep_or_shallow"],
        8: forms["slope_sqrt_three"],
        10: forms["slope_one_over_sqrt_three"],
        11: forms["negative_steep_or_shallow"],
    }
    rational_targets = {0, 3, 6, 9}
    assert set(irrational_target_forms) | rational_targets == set(range(12))
    assert set(irrational_target_forms).isdisjoint(rational_targets)
    return len(forms)


def integer_partitions(total, largest=None):
    if total == 0:
        yield ()
        return
    if largest is None or largest > total:
        largest = total
    for first in range(largest, 0, -1):
        for rest in integer_partitions(total - first, first):
            yield (first,) + rest


def class_labels(parts):
    labels = []
    for class_index, size in enumerate(parts):
        labels.extend([class_index] * size)
    return labels


def matching_edges(k):
    if k % 2:
        step = (k - 1) // 2
        edges = {tuple(sorted((i, (i + step) % k))) for i in range(k)}
        weight = Fraction(1, 2)
    else:
        step = k // 2
        edges = {(i, i + step) for i in range(step)}
        weight = Fraction(1)
    return edges, weight


def check_class_matchings(max_k=24):
    partition_count = 0
    for k in range(2, max_k + 1):
        for parts in integer_partitions(k):
            if max(parts) > k // 2:
                continue
            partition_count += 1
            labels = class_labels(parts)
            edges, weight = matching_edges(k)
            assert all(labels[a] != labels[b] for a, b in edges)

            loads = [Fraction(0) for _ in range(k)]
            for a, b in edges:
                loads[a] += weight
                loads[b] += weight
            assert loads == [Fraction(1)] * k
            assert sum((weight for _ in edges), Fraction(0)) == Fraction(k, 2)
    return partition_count


def main():
    pair_count = check_paley_core()
    form_count = check_quadratic_forms()
    partition_count = check_class_matchings()
    print(
        "verified",
        pair_count,
        "Paley pair isolators,",
        form_count,
        "quadratic derivative bounds, and",
        partition_count,
        "class-partition matchings",
    )


if __name__ == "__main__":
    main()
