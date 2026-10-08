"""Exact checks for partition_kernel_phase_groups.md.

Run with python3 -B docs/check_partition_kernel_phase_groups.py.
The all-b proof and the Diophantine separation are in the prose note.
This checks the full row spaces, not merely a selected-column submatrix.
"""

from fractions import Fraction

from check_rank_one_kernel_relative_roth import paley_core_11


def rank(matrix):
    work = [[Fraction(x) for x in row] for row in matrix]
    pivot_row = 0
    for col in range(len(work[0])):
        pivot = next(
            (i for i in range(pivot_row, len(work)) if work[i][col]), None
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [x / value for x in work[pivot_row]]
        for i in range(pivot_row + 1, len(work)):
            value = work[i][col]
            if value:
                work[i] = [
                    x - value * y for x, y in zip(work[i], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def walk_matrix(groups):
    core = paley_core_11()
    state = [0] * groups
    states = [state[:]]
    increments = {}
    for h in range(1, 11):
        for g in range(groups):
            state[g] = h
            states.append(state[:])
            increments[g, h] = len(states) - 1
    rows = [sum((core[a] for a in state), []) for state in states]
    return core, rows, increments


def check_walk(groups):
    core, rows, increments = walk_matrix(groups)
    row_count = 10 * groups + 1
    col_count = 11 * groups
    assert len(rows) == row_count
    assert all(len(row) == col_count for row in rows)

    # Consecutive differences and differences from row zero have equal span.
    differences = [
        [a - b for a, b in zip(rows[t], rows[t - 1])]
        for t in range(1, row_count)
    ]
    assert rank(differences) == 10 * groups
    for g in range(groups):
        indicator = [int(j // 11 == g) for j in range(col_count)]
        assert all(sum(a * b for a, b in zip(row, indicator)) == 0
                   for row in differences)
        assert all(sum(row[11 * g:11 * (g + 1)]) == -1 for row in rows)

    # Verify every certificate against all 11b columns, including other groups.
    certificate_count = 0
    maximum_load = Fraction(0)
    for g in range(groups):
        for a in range(11):
            for c in range(a):
                local = [Fraction(core[i][a] - core[i][c], 12)
                         for i in range(11)]
                assert sum(local) == 0
                assert sum(abs(x) for x in local) == 1
                coefficients = [Fraction(0)] * row_count
                for h in range(1, 11):
                    tail = sum(local[h:])
                    after = increments[g, h]
                    coefficients[after] += tail
                    coefficients[after - 1] -= tail
                assert sum(coefficients) == 0
                assert all((6 * x).denominator == 1 for x in coefficients)
                load = sum(abs(x) for x in coefficients)
                assert load <= 20
                maximum_load = max(maximum_load, load)
                expected = [Fraction(int(j == 11 * g + a)
                                     - int(j == 11 * g + c))
                            for j in range(col_count)]
                actual = [sum(coefficients[i] * rows[i][j]
                              for i in range(row_count))
                          for j in range(col_count)]
                assert actual == expected
                certificate_count += 1

    # Matching weights and exponent conversion are exact rational identities.
    weight = Fraction(11 * groups, 2)
    delta_exponent = Fraction(2, 1) / weight
    arc_exponent = 1 - delta_exponent
    assert arc_exponent == 1 - Fraction(4, 11 * groups)
    endpoint_exponent = 1 / (arc_exponent - Fraction(1, 2))
    assert endpoint_exponent == Fraction(22 * groups, 11 * groups - 8)
    lower_radius_exponent = Fraction(11 * groups, 2)
    assert lower_radius_exponent * (arc_exponent - Fraction(1, 2)) == (
        Fraction(11 * groups - 8, 4)
    )
    return certificate_count, maximum_load


def main():
    total = 0
    for groups in (1, 2, 3, 4):
        count, load = check_walk(groups)
        total += count
        print("b={}: kernel dimension {}, {} certificates, max l1={}".format(
            groups, groups, count, load
        ))
    print("PASS: {} full-space certificates and all stated exponent identities.".format(total))


if __name__ == "__main__":
    main()
