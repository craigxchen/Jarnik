"""Exact one-section boundary slopes on the five-row del Pezzo surface."""

from fractions import Fraction
from itertools import combinations

from check_five_row_boundary_node_height import linear_coefficients
from check_five_row_del_pezzo_arithmetic import graphs
from check_smaller_cut_relation_rank import restrict_graph


PAIRS = tuple(combinations(range(5), 2))
GRAPH_INDICES = (0, 1, 2, 3, 4, 5)

# The expected (alpha,beta) rows, in the six graph coordinates G0,...,G5.
EXPECTED = {
    (0, 1): ((0, 0, 0, -1, -1, 0), (0, 0, 0, 1, 0, -1)),
    (0, 2): ((0, -1, -1, 0, 0, 0), (0, 1, 0, 0, 0, -1)),
    (0, 3): ((-1, 0, -1, 0, 0, 0), (1, 0, 0, 0, -1, 0)),
    (0, 4): ((-1, -1, 0, 0, 0, 0), (1, 0, 0, -1, 0, 0)),
    (1, 2): ((0, 1, 1, 1, 1, 0), (0, -1, 0, -1, 0, 0)),
    (1, 3): ((1, 0, 1, -1, 0, 1), (-1, 0, 0, 1, 0, 0)),
    (1, 4): ((1, 1, 0, 0, -1, -1), (-1, 0, 0, 0, 1, 0)),
    (2, 3): ((-1, -1, 0, 0, 1, 1), (1, 1, 0, 0, 0, 0)),
    (2, 4): ((-1, 0, -1, 1, 0, -1), (1, 0, 1, 0, 0, 0)),
    (3, 4): ((0, -1, -1, -1, -1, 0), (0, 1, 1, 0, 0, 0)),
}


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


def coefficient_vector(graph, inside):
    return linear_coefficients(restrict_graph(graph, set(inside), 5))


def node_form(outside, remaining):
    """Coefficient of z_remaining as a form in (alpha,beta)."""
    position = outside.index(remaining)
    return ((1, 0), (0, 1), (-1, -1))[position]


# These forms reproduce the displayed identities in the note.  Each pair is
# (coefficient of alpha, coefficient of beta) on the corresponding line.
NODE_RELATIONS = (
    (((0, 1), (1, 1)), ((2, 3), (1, 1))),
    (((0, 1), (0, 1)), ((2, 4), (-1, -1))),
    (((0, 1), (1, 0)), ((3, 4), (-1, -1))),
    (((0, 2), (1, 1)), ((1, 3), (1, 1))),
    (((0, 2), (0, 1)), ((1, 4), (-1, -1))),
    (((0, 2), (1, 0)), ((3, 4), (0, 1))),
    (((0, 3), (1, 1)), ((1, 2), (1, 1))),
    (((0, 3), (0, 1)), ((1, 4), (0, 1))),
    (((0, 3), (1, 0)), ((2, 4), (0, 1))),
    (((0, 4), (-1, -1)), ((1, 2), (0, 1))),
    (((0, 4), (0, 1)), ((1, 3), (0, 1))),
    (((0, 4), (1, 0)), ((2, 3), (0, 1))),
    (((1, 2), (1, 0)), ((3, 4), (1, 0))),
    (((1, 3), (1, 0)), ((2, 4), (1, 0))),
    (((1, 4), (1, 0)), ((2, 3), (1, 0))),
)


def main():
    all_graphs = list(graphs())
    basis = [all_graphs[j] for j in GRAPH_INDICES]
    # Verify independence using all polynomial coefficients at the empty cut.
    monomials = sorted(set().union(*(restrict_graph(g, set(), 5) for g in basis)))
    assert rank([[restrict_graph(g, set(), 5).get(m, 0) for m in monomials]
                for g in basis]) == 6

    actual = {}
    for inside in PAIRS:
        outside = tuple(j for j in range(5) if j not in inside)
        alpha, beta = [], []
        for graph in basis:
            row = coefficient_vector(graph, inside)
            alpha.append(row[outside[0]])
            beta.append(row[outside[1]])
        actual[inside] = (tuple(alpha), tuple(beta), outside)
        assert (tuple(alpha), tuple(beta)) == EXPECTED[inside]
        assert all(a + b + row[outside[2]] == 0
                   for a, b, row in zip(alpha, beta,
                                        [coefficient_vector(g, inside) for g in basis]))

    # The shared-chart node coefficient is opposite on the two incident lines.
    node_count = 0
    for first, second in combinations(PAIRS, 2):
        if set(first).isdisjoint(second):
            remaining = next(j for j in range(5) if j not in first + second)
            for coordinate in range(6):
                a, b, outside = actual[first]
                c, d, outside2 = actual[second]
                fa, fb = node_form(outside, remaining)
                fc, fd = node_form(outside2, remaining)
                assert fa * a[coordinate] + fb * b[coordinate] == -(
                    fc * c[coordinate] + fd * d[coordinate]
                )
            node_count += 1
    assert node_count == 15
    # Also test the exact signs of the identities printed in the note,
    # rather than only the equivalent opposite-coefficient statement.
    for (first, first_form), (second, second_form) in NODE_RELATIONS:
        a, b, _ = actual[first]
        c, d, _ = actual[second]
        fa, fb = first_form
        fc, fd = second_form
        assert all(
            fa * a[j] + fb * b[j] + fc * c[j] + fd * d[j] == 0
            for j in range(6)
        )
        # The relation just checked is the unique linear relation: the
        # four forms span a rank-three hyperplane in q-space.
        assert rank([a, b, c, d]) == 3
    assert len(NODE_RELATIONS) == 15
    print("PASS: six graph coordinates span; all ten exact (alpha,beta) rows match.")
    print("PASS: all 15 disjoint boundary pairs satisfy the shared-node coefficient identity.")
    print("PASS: each incident four-form tuple has exact rank three, so the node relation is unique.")
    print("No second section, numerical zero, or scale-free slope-ratio theorem is asserted.")


if __name__ == "__main__":
    main()
