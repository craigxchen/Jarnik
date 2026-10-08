"""Finite identity audit for the weighted Hadamard-core replacement lemma.

The theorem itself is proved for every Hadamard order in the companion note.
This checker exhausts replacement columns at Sylvester orders 4 and 8 and
verifies a weight-independent obstruction for every non-core replacement.
"""

from itertools import product


def sylvester(order):
    h = [[1]]
    while len(h) < order:
        h = [row + row for row in h] + [row + [-x for x in row] for row in h]
    return h


def audit(order):
    h = sylvester(order)
    columns = [[h[i][j] for i in range(order)] for j in range(order)]
    assert columns[0] == [1] * order
    for omitted in range(1, order):
        chi = columns[omitted]
        a = [i for i in range(order) if chi[i] == 1]
        b = [i for i in range(order) if chi[i] == -1]
        old = [columns[j] for j in range(1, order) if j != omitted]
        for col in old:
            assert sum(col[i] for i in a) == 0
            assert sum(col[i] for i in b) == 0
        for f in product((1, -1), repeat=order):
            fa = sum(f[i] for i in a)
            fb = sum(f[i] for i in b)
            # For each cross-part row, summing its Gram entries cancels
            # every retained column for arbitrary positive old weights.
            cross_row_sums = [f[i] * fb for i in a] + [f[i] * fa for i in b]
            if tuple(f) in (tuple(chi), tuple(-x for x in chi)):
                assert all(z < 0 for z in cross_row_sums)
            else:
                assert any(z >= 0 for z in cross_row_sums)
    return (order - 1) * 2 ** order


if __name__ == "__main__":
    for order in (4, 8):
        print("order", order, "replacement cases", audit(order), "passed")
