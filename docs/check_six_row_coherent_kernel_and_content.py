"""Symbolic certificate for the complete six-row coherent cut kernel.

All polynomial identities use six independent affine variables, then extend
to binary rows by separate homogeneity. No sampled-rank inference is used.
"""

from itertools import combinations, permutations
from math import comb

from check_near_balanced_invariant_congruences import ZERO, add, mul, scale, variable


def product(polynomials):
    answer = {ZERO: 1}
    for polynomial in polynomials:
        answer = mul(answer, polynomial)
    return answer


def delta(i, j):
    return add(variable(j), scale(variable(i), -1))


MATCHINGS = (
    ((0, 1), (2, 3), (4, 5)),
    ((0, 1), (2, 5), (3, 4)),
    ((0, 3), (1, 2), (4, 5)),
    ((0, 5), (1, 2), (3, 4)),
    ((0, 5), (1, 4), (2, 3)),
)
PARTITIONS = [(J, tuple(i for i in range(6) if i not in J))
              for J in combinations(range(6), 3) if 0 in J]
CUTS = list(combinations(range(6), 2))
TREE = ((0, 1), (3, 4), (1, 2), (2, 3))
F_EDGES = ((0, 1), (4, 5), (2, 5), (0, 5), (3, 4))
UV_ROWS = (
    ((1,1,0,0,0), (-1,0,0,0,0)),
    ((-1,-1,-1,-1,0), (1,0,1,0,0)),
    ((1,1,1,1,1), (0,0,-1,0,0)),
    ((-1,0,-1,-1,-1), (0,0,1,1,0)),
    ((0,0,0,1,1), (0,0,0,-1,0)),
    ((0,0,1,1,0), (0,0,-1,0,0)),
    ((-1,0,-1,-1,-1), (1,0,1,0,0)),
    ((1,1,1,1,1), (-1,-1,-1,-1,0)),
    ((-1,-1,0,-1,-1), (0,1,0,1,0)),
    ((1,0,0,0,1), (-1,0,0,0,0)),
    ((-1,-1,0,-1,-1), (1,1,0,0,0)),
    ((1,1,1,1,1), (0,-1,0,0,0)),
    ((0,1,0,1,0), (0,-1,0,0,0)),
    ((-1,-1,-1,-1,0), (1,1,0,0,0)),
    ((1,0,1,0,0), (-1,0,0,0,0)),
)

# Columns express G_A,...,G_E in the ten triangle products.
CONVERSION = (
    (2, 1, 1, 1, 1),
    (-1, -1, -1, -1, 0),
    (1, 0, 1, 0, 0),
    (0, 0, 0, 0, 0),
    (1, 1, 0, 0, 0),
    (-1, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0),
)


def triangle_cut(partition, S):
    factors = []
    for a, b, c in partition:
        inside = set((a, b, c)) & set(S)
        if len(inside) != 1:
            return {}
        i = next(iter(inside))
        factors.append({a: delta(b, c),
                        b: scale(delta(a, c), -1),
                        c: delta(a, b)}[i])
    return product(factors)


def coherent_matrix():
    matrix = []
    for S in CUTS:
        triangles = [triangle_cut(partition, S) for partition in PARTITIONS]
        row = []
        for column in range(5):
            value = {}
            for i in range(10):
                value = add(value, scale(triangles[i], CONVERSION[i][column]))
            row.append(value)
        matrix.append(row)
    return matrix


def determinant(matrix):
    n = len(matrix)
    answer = {}
    for ordering in permutations(range(n)):
        inversions = sum(ordering[i] > ordering[j]
                         for i in range(n) for j in range(i+1, n))
        term = product(matrix[i][ordering[i]] for i in range(n))
        answer = add(answer, scale(term, (-1)**inversions))
    return answer


def symbolic_checks():
    M = [product(delta(i, j) for i, j in matching) for matching in MATCHINGS]
    matrix = coherent_matrix()
    for S, row, (u, v) in zip(CUTS, matrix, UV_ROWS):
        a, b, c, d = [i for i in range(6) if i not in S]
        U = mul(delta(a, b), delta(c, d))
        V = mul(delta(a, c), delta(b, d))
        assert row == [add(scale(U, left), scale(V, right))
                       for left, right in zip(u, v)]
    for row in matrix:
        answer = {}
        for entry, matching in zip(row, M):
            answer = add(answer, mul(entry, matching))
        assert not answer

    # Verify the four sparse rows without polynomial division.
    sparse = (
        ((0, 1), 0, 1),
        ((3, 4), 1, 3),
        ((1, 2), 2, 3),
        ((2, 3), 0, 4),
    )
    for S, left, right in sparse:
        row = matrix[CUTS.index(S)]
        for column in range(5):
            expected = scale(M[right], -1) if column == left else M[left] if column == right else {}
            assert mul(delta(*S), row[column]) == expected

    selected = [matrix[CUTS.index(S)] for S in TREE]
    F = product(delta(i, j) for i, j in F_EDGES)
    for omitted in range(5):
        minor = determinant([[entry for j, entry in enumerate(row) if j != omitted]
                             for row in selected])
        assert scale(minor, (-1)**omitted) == mul(F, M[omitted])
    print("PASS: all15 symbolic kernel equations, four sparse rows, and all five tree cofactors.")


def content_checks():
    active = []
    total = 0
    tree = [set(S) for S in TREE]
    for size in range(7):
        for labels in combinations(range(6), size):
            T = set(labels)
            match_order = min(sum(i in T and j in T for i, j in matching)
                              for matching in MATCHINGS)
            assert match_order == max(0, size-3)
            F_order = sum(i in T and j in T for i, j in F_EDGES)
            row_order = sum(max(0, len(T-S)-2) for S in tree)
            exponent = F_order + match_order - row_order
            assert exponent in (0, 1)
            if exponent:
                active.append("".join(str(i+1) for i in labels))
            total += exponent
    assert active == ["12", "16", "36", "45", "56", "124", "136", "245", "356",
                      "1234", "1236", "1245", "2345", "3456"]
    assert total == 14
    assert sum(comb(6, a)*max(0, a-3) for a in range(7)) == 30
    for S in tree:
        content = sum(max(0, len(set(T)-S)-2)
                      for a in range(7) for T in combinations(range(6), a))
        assert content == 24
    assert 3*16-30 == 18
    assert 4*(2*16-24)-14 == 18
    print("PASS: all64 exact cut orders; selected normalized cofactor content is14w,")
    print("      primitive kernel height is18w, above the8w actual-cut exclusion threshold.")


if __name__ == "__main__":
    symbolic_checks()
    content_checks()
