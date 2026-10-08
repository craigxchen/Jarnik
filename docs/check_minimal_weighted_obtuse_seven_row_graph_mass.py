"""Exact finite audit of the (2,1^6) seven-row graph obstruction."""

from itertools import combinations


def audit():
    vertices = tuple(range(6))
    all_edges = tuple(combinations(vertices, 2))
    low_degree = 0
    two_regular = 0
    for chosen in combinations(all_edges, 6):
        degree = [sum(v in e for e in chosen) for v in vertices]
        if min(degree) <= 1:
            low_degree += 1
            continue
        assert degree == [2] * 6
        two_regular += 1
        # For each selected endpoint pair, the edge itself contributes
        # zero to D(i,k); every adjacent edge contributes its full weight.
        # Each edge is adjacent to exactly two other selected edges.
        for edge in chosen:
            assert sum(len(set(edge) & set(other)) == 1
                       for other in chosen if other != edge) == 2
    assert low_degree + two_regular == 5005
    return low_degree, two_regular


if __name__ == "__main__":
    print("low-degree or two-regular six-edge graphs:", audit())
