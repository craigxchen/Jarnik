"""Exact Petersen circulation and gluing-rank checks for five rows."""

from fractions import Fraction
from itertools import combinations
from math import gcd

from check_five_row_boundary_node_height import linear_coefficients
from check_five_row_del_pezzo_arithmetic import graphs
from check_smaller_cut_relation_rank import restrict_graph


PAIRS = tuple(combinations(range(5), 2))
EDGES = tuple(
    (i, j) for i, S in enumerate(PAIRS) for j, T in enumerate(PAIRS)
    if i < j and set(S).isdisjoint(T)
)
SIX_GRAPH_BASIS = tuple(list(graphs())[:6])
ALL_GRAPHS = tuple(graphs())


def rank(rows):
    rows = [[Fraction(value) for value in row] for row in rows]
    pivots = {}
    for row in rows:
        for pivot, previous in sorted(pivots.items()):
            if row[pivot]:
                factor = row[pivot]
                row = [a - factor * b for a, b in zip(row, previous)]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            factor = row[pivot]
            pivots[pivot] = [value / factor for value in row]
    return len(pivots)


def determinant(matrix):
    matrix = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(matrix)):
        pivot = next((j for j in range(column, len(matrix))
                      if matrix[j][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            matrix[pivot], matrix[column] = matrix[column], matrix[pivot]
            result = -result
        value = matrix[column][column]
        result *= value
        for row in range(column + 1, len(matrix)):
            factor = matrix[row][column] / value
            for entry in range(column + 1, len(matrix)):
                matrix[row][entry] -= factor * matrix[column][entry]
    return result


def coefficient_table(graph_basis):
    table = {}
    for pair in PAIRS:
        table[pair] = [
            linear_coefficients(restrict_graph(graph, set(pair), 5))
            for graph in graph_basis
        ]
    return table


def main():
    assert len(PAIRS) == 10 and len(EDGES) == 15
    table = coefficient_table(ALL_GRAPHS)

    incidence = [[0] * len(EDGES) for _ in PAIRS]
    flow_matrix = []
    for edge_index, (first_index, second_index) in enumerate(EDGES):
        first, second = PAIRS[first_index], PAIRS[second_index]
        remaining = next(i for i in range(5) if i not in first + second)
        outside = tuple(i for i in range(5) if i not in first)
        outside_second = tuple(i for i in range(5) if i not in second)
        assert remaining in outside and remaining in outside_second
        flow_column = []
        for coordinate in range(len(ALL_GRAPHS)):
            left = table[first][coordinate][remaining]
            right = table[second][coordinate][remaining]
            assert left == -right
            flow_column.append(left)
        flow_matrix.append(flow_column)
        incidence[first_index][edge_index] = 1
        incidence[second_index][edge_index] = -1

    assert rank(incidence) == 9
    assert rank(flow_matrix) == 6
    assert all(
        sum(incidence[vertex][edge] * flow_matrix[edge][coordinate]
            for edge in range(15)) == 0
        for vertex in range(10) for coordinate in range(len(ALL_GRAPHS))
    )
    # Each graph generator is a primitive signed cycle flow: the ten
    # triangle-times-squared-edge graphs give 6-cycles, and the twelve
    # pentagons give 5-cycles.  This also records the visible parity split.
    for coordinate in range(len(ALL_GRAPHS)):
        support = [edge for edge in range(15) if flow_matrix[edge][coordinate]]
        expected_length = 6 if coordinate < 10 else 5
        assert len(support) == expected_length
        assert all(abs(flow_matrix[edge][coordinate]) == 1 for edge in support)
        adjacency = {vertex: [] for vertex in range(10)}
        for edge in support:
            first, second = EDGES[edge]
            adjacency[first].append(second)
            adjacency[second].append(first)
        assert all(len(neighbors) in (0, 2) for neighbors in adjacency.values())
        assert sum(bool(neighbors) for neighbors in adjacency.values()) == expected_length
        reached = {support and EDGES[support[0]][0]}
        while True:
            expanded = reached | {
                neighbor for vertex in reached for neighbor in adjacency[vertex]
            }
            if expanded == reached:
                break
            reached = expanded
        assert len(reached) == expected_length
        assert sum(flow_matrix[edge][coordinate] for edge in support) % 2 == expected_length % 2

    # The first six graph columns form the old index-two sublattice.  The
    # complete 22-column graph map contains an actual unimodular flow minor.
    first_six = [column[:6] for column in flow_matrix]
    first_six_transpose = [[first_six[row][coordinate] for coordinate in range(6)]
                           for row in range(15)]
    first_six_minors = [
        determinant([first_six_transpose[row] for row in rows])
        for rows in combinations(range(15), 6)
    ]
    index = 0
    for minor in first_six_minors:
        index = gcd(index, abs(minor.numerator))
    assert index == 2 and sum(minor != 0 for minor in first_six_minors) == 2000

    witness_rows = (0, 7, 8, 9, 13, 14)
    witness_columns = (0, 1, 2, 3, 4, 10)
    witness = determinant([
        [flow_matrix[row][column] for column in witness_columns]
        for row in witness_rows
    ])
    assert witness == -1

    # Build the 15 gluing equations in the 20 alpha/beta coordinates.
    gluing = []
    for first_index, second_index in EDGES:
        first, second = PAIRS[first_index], PAIRS[second_index]
        remaining = next(i for i in range(5) if i not in first + second)
        row = [0] * 20
        for index, pair in ((first_index, first), (second_index, second)):
            outside = tuple(i for i in range(5) if i not in pair)
            position = outside.index(remaining)
            if position == 0:
                row[2 * index] += 1
            elif position == 1:
                row[2 * index + 1] += 1
            else:
                row[2 * index] -= 1
                row[2 * index + 1] -= 1
        gluing.append(row)
    assert rank(gluing) == 14

    # Coefficientwise global identity in common z_1,...,z_5 variables.
    for coordinate in range(6):
        assert all(
            sum(table[pair][coordinate][variable] for pair in PAIRS) == 0
            for variable in range(5)
        )

    q = (1, 2, 3, 4, 5, 6) + (0,) * 16
    common_forms = [
        [sum(q[coordinate] * table[pair][coordinate][variable]
             for coordinate in range(6)) for variable in range(5)]
        for pair in PAIRS
    ]
    assert rank(common_forms) == 4
    assert all(sum(form) == 0 for form in common_forms)

    print("PASS: Petersen incidence rank 9; q-to-flow rank 6; circulation exact.")
    print("PASS: all 10 triangle-square flows are signed 6-cycles; all 12 pentagon flows are signed 5-cycles.")
    print("PASS: first-six flow sublattice has index 2 (2000 nonzero maximal minors).")
    print("PASS: all-22 graph flow image contains a unimodular minor, so its index is 1.")
    print("PASS: 15 gluing equations have rank 14; sum_S F_S(z)=0 coefficientwise.")
    print("PASS: fixed q=(1,...,6) gives rank-4 common-variable forms.")
    print("No Gaussian-modulus multiplication or uniform arc conclusion is asserted.")


if __name__ == "__main__":
    main()
