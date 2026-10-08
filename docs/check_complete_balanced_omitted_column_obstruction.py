"""Finite checks for complete_balanced_omitted_column_obstruction.md.

Standard library only. This checks exact rational row-space membership and
the balanced-column counting formulas. It does not model Gaussian phases or
claim an endpoint realization.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb


def rank(matrix):
    work = [[Fraction(x) for x in row] for row in matrix]
    if not work:
        return 0
    pivot = 0
    for col in range(len(work[0])):
        found = next((i for i in range(pivot, len(work)) if work[i][col]), None)
        if found is None:
            continue
        work[pivot], work[found] = work[found], work[pivot]
        scale = work[pivot][col]
        work[pivot] = [x / scale for x in work[pivot]]
        for i in range(len(work)):
            if i != pivot and work[i][col]:
                scale = work[i][col]
                work[i] = [x - scale * y for x, y in zip(work[i], work[pivot])]
        pivot += 1
        if pivot == len(work):
            break
    return pivot


def nullspace(matrix):
    """A rational basis for the right kernel of a row matrix."""
    work = [[Fraction(x) for x in row] for row in matrix]
    columns = len(work[0])
    pivot_row = 0
    pivot_columns = []
    for col in range(columns):
        found = next((i for i in range(pivot_row, len(work))
                      if work[i][col]), None)
        if found is None:
            continue
        work[pivot_row], work[found] = work[found], work[pivot_row]
        scale = work[pivot_row][col]
        work[pivot_row] = [x / scale for x in work[pivot_row]]
        for i in range(len(work)):
            if i != pivot_row and work[i][col]:
                scale = work[i][col]
                work[i] = [x - scale * y for x, y in zip(work[i], work[pivot_row])]
        pivot_columns.append(col)
        pivot_row += 1
    free_columns = [j for j in range(columns) if j not in pivot_columns]
    basis = []
    for free in free_columns:
        vector = [Fraction(0)] * columns
        vector[free] = 1
        for row, pivot in reversed(list(enumerate(pivot_columns))):
            vector[pivot] = -work[row][free]
        basis.append(vector)
    return basis


def in_rowspace(kernel_basis, vector):
    return all(sum(a * b for a, b in zip(q, vector)) == 0
               for q in kernel_basis)


def kernel_projection_rank(kernel_basis, indices):
    return rank([[q[j] for j in indices] for q in kernel_basis])


def balanced_sign_matrix(m):
    k = m // 2
    columns = tuple(combinations(range(m), k))
    rows = [[1 if i in column else -1 for column in columns] for i in range(m)]
    differences = [[x - y for x, y in zip(rows[i], rows[0])]
                   for i in range(1, m)]
    return columns, rows, differences


def variable_count(m, t):
    k = m // 2
    b = comb(m, k)
    if t <= k:
        return b - 2 * comb(m - t, k - t)
    return b


def reduced_patterns(m, t):
    """Restricted nonconstant cuts, modulo complement."""
    k = m // 2
    patterns = set()
    for signs in product((-1, 1), repeat=t):
        h = sum(x == 1 for x in signs)
        if 0 < h < t and max(0, t - k) <= h <= min(t, k):
            patterns.add(min(tuple(signs), tuple(-x for x in signs)))
    return sorted(patterns)


def reduced_count(m, t):
    k = m // 2
    lo, hi = max(1, t - k), min(t - 1, k)
    return sum(comb(t, h) for h in range(lo, hi + 1)) // 2


def check_reduced_patterns(m):
    total = 0
    full_columns, _, _ = balanced_sign_matrix(m)
    for t in range(1, m + 1):
        patterns = reduced_patterns(m, t)
        # Independently restrict the actual original balanced columns, rather
        # than only comparing two uses of the feasible-size formula.
        actual_patterns = set()
        for column in full_columns:
            signs = tuple(1 if i in column else -1 for i in range(t))
            if len(set(signs)) > 1:
                actual_patterns.add(min(signs, tuple(-x for x in signs)))
        assert patterns == sorted(actual_patterns)
        assert len(patterns) == reduced_count(m, t)
        total += len(patterns)
        if not patterns:
            continue
        rows = [list(column) for column in zip(*patterns)]
        differences = [[x - y for x, y in zip(rows[i], rows[0])]
                       for i in range(1, t)]
        assert rank(differences) == min(t - 1, len(patterns))
        if t >= 4:
            kernel_basis = nullspace(differences)
            # No nonzero vector with arbitrary coefficients supported on one
            # or two reduced classes.
            for j in range(len(patterns)):
                assert kernel_projection_rank(kernel_basis, (j,)) == 1
                for h in range(j):
                    assert kernel_projection_rank(kernel_basis, (h, j)) == 2
        if t >= 4:
            assert len(patterns) > t
    return total


def check_matrix(m):
    columns, rows, differences = balanced_sign_matrix(m)
    b = len(columns)
    assert rank(differences) == m - 1
    kernel_basis = nullspace(differences)
    assert all(sum(row[j] for row in rows) == 0 for j in range(b))
    assert rank(rows) == rank(differences) == m - 1

    # Every nonzero vector with support at most two misses the row space.
    tested = 0
    for j in range(b):
        assert kernel_projection_rank(kernel_basis, (j,)) == 1
        tested += 1
        for h in range(j):
            assert kernel_projection_rank(kernel_basis, (h, j)) == 2
            tested += 1

    # Count original varying columns before the regrouping operation.
    for t in range(1, m + 1):
        variable = sum(
            len({i in column for i in range(t)}) == 2 for column in columns
        )
        assert variable == variable_count(m, t)
        if t >= 2:
            assert variable > t
    return b, tested


def main():
    total_columns = total_vectors = total_reduced = 0
    for m in (6, 8):
        columns, tested = check_matrix(m)
        total_columns += columns
        total_vectors += tested
        total_reduced += check_reduced_patterns(m)
        print("M={}: {} balanced columns, {} support tests".format(
            m, columns, tested
        ))
    print("PASS: {} columns, {} original support tests, and {} reduced classes.".format(
        total_columns, total_vectors, total_reduced
    ))
    print("No endpoint realization or uniform counting theorem is asserted.")


if __name__ == "__main__":
    main()
