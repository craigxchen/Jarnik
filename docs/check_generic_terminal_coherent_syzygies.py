"""Exact finite certificates for generic terminal coherent syzygies.

All ranks are computed by modular row reduction over F_101.  A nonzero
minor modulo 101 is a nonzero integer minor, so the reported lower ranks
hold over Q.  This checker does not assert constancy at every configuration.
"""

from itertools import combinations
import sys

import numpy as np


PRIME = 101
LABELS = tuple(range(12))
SIX_SETS = tuple(combinations(LABELS, 6))
CUTS = tuple(combinations(LABELS, 4))
MATCHINGS = (
    ((0, 1), (2, 3), (4, 5)),
    ((0, 1), (2, 5), (3, 4)),
    ((0, 3), (1, 2), (4, 5)),
    ((0, 5), (1, 2), (3, 4)),
    ((0, 5), (1, 4), (2, 3)),
)


def modular_rank(values):
    """Integer row rank modulo PRIME; int64 products here are < 2**63."""
    matrix = np.array(values, dtype=np.int64, copy=True) % PRIME
    nrows, ncols = matrix.shape
    rank = 0
    for column in range(ncols):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        selected = rank + int(candidates[0])
        matrix[[rank, selected]] = matrix[[selected, rank]]
        inverse = pow(int(matrix[rank, column]), -1, PRIME)
        matrix[rank, column:] = matrix[rank, column:] * inverse % PRIME
        factors = matrix[rank + 1:, column].copy()
        matrix[rank + 1:, column:] = (
            matrix[rank + 1:, column:]
            - factors[:, None] * matrix[rank, column:]
        ) % PRIME
        rank += 1
        if rank == nrows:
            break
    return rank


def standard_tableaux(rows=((), (), ())):
    """Increasing 3-by-4 tableaux, used only to supply 462 K_2 vectors."""
    next_label = sum(map(len, rows))
    if next_label == 12:
        yield rows
        return
    for row in range(3):
        if len(rows[row]) == 4 or (row and len(rows[row - 1]) <= len(rows[row])):
            continue
        extended = list(rows)
        extended[row] += (next_label,)
        yield from standard_tableaux(tuple(extended))


TABLEAUX = tuple(standard_tableaux())
TRIANGLE_COLUMNS = tuple(tuple(zip(*tableau)) for tableau in TABLEAUX)


def coherent_matrix(directions):
    """Rows indexed by four-label cuts, columns by four-triangle graphs."""
    result = np.zeros((len(CUTS), len(TRIANGLE_COLUMNS)), dtype=np.int64)
    for cut_index, cut in enumerate(CUTS):
        inside = set(cut)
        for column_index, triangles in enumerate(TRIANGLE_COLUMNS):
            value = 1
            for triangle in triangles:
                hits = [j for j, label in enumerate(triangle) if label in inside]
                if len(hits) != 1:
                    value = 0
                    break
                outside_pair = [label for label in triangle if label not in inside]
                orientation = -1 if hits[0] == 1 else 1
                value *= orientation * (
                    directions[outside_pair[1]] - directions[outside_pair[0]]
                )
            result[cut_index, column_index] = value % PRIME
    return result


def test_vectors():
    """462 deterministic ternary configurations; no random library state."""
    state = 20261001
    entries = []
    for _ in range(462 * 12 * 3):
        state = (1664525 * state + 1013904223) % (2**32)
        entries.append(state % PRIME)
    return np.array(entries, dtype=np.int64).reshape(462, 12, 3)


def triple_determinants(vectors):
    answer = {}
    for i, j, k in combinations(LABELS, 3):
        a, b, c = vectors[:, i, :], vectors[:, j, :], vectors[:, k, :]
        answer[(i, j, k)] = (
            a[:, 0] * (b[:, 1] * c[:, 2] - b[:, 2] * c[:, 1])
            - a[:, 1] * (b[:, 0] * c[:, 2] - b[:, 2] * c[:, 0])
            + a[:, 2] * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
        ) % PRIME
    return answer


def gradient_values(six, determinants):
    """The exact T_123,...,T_135 formulas for G_A,...,G_E."""
    def pair_product(triple):
        triple = set(triple)
        other = tuple(i for i in range(6) if i not in triple)
        left = tuple(six[i] for i in sorted(triple))
        right = tuple(six[i] for i in other)
        return determinants[left] * determinants[right] % PRIME

    t123 = pair_product((0, 1, 2))
    t124 = pair_product((0, 1, 3))
    t125 = pair_product((0, 1, 4))
    t134 = pair_product((0, 2, 3))
    t135 = pair_product((0, 2, 4))
    return np.stack((
        (2*t123 - t124 + t125 + t134 - t135) % PRIME,
        (t123 - t124 + t134) % PRIME,
        (t123 - t124 + t125) % PRIME,
        (t123 - t124) % PRIME,
        t123,
    ), axis=1)


def matching_values(six, directions):
    def bracket(i, j):
        return directions[six[j]] - directions[six[i]]
    return np.array([
        np.prod([bracket(i, j) for i, j in matching]) % PRIME
        for matching in MATCHINGS
    ], dtype=np.int64)


def check_local_six_row_incidence(directions):
    """Check all 924 local polar vectors against their 15 coherent cuts."""
    selected_triples = ((0, 1, 2), (0, 1, 3), (0, 1, 4),
                        (0, 2, 3), (0, 2, 4))
    for six in SIX_SETS:
        matching = matching_values(six, directions)
        for local_cut in combinations(range(6), 2):
            inside = set(local_cut)
            triangle_pairs = []
            for triple in selected_triples:
                complement = tuple(i for i in range(6) if i not in triple)
                factors = []
                for triangle in (triple, complement):
                    hits = [j for j, label in enumerate(triangle) if label in inside]
                    if len(hits) != 1:
                        factors = [0]
                        break
                    outside = [label for label in triangle if label not in inside]
                    orientation = -1 if hits[0] == 1 else 1
                    factors.append(orientation * (
                        directions[six[outside[1]]] - directions[six[outside[0]]]
                    ))
                triangle_pairs.append(np.prod(factors) % PRIME)
            t123, t124, t125, t134, t135 = triangle_pairs
            gradient = np.array((
                2*t123 - t124 + t125 + t134 - t135,
                t123 - t124 + t134,
                t123 - t124 + t125,
                t123 - t124,
                t123,
            ), dtype=np.int64) % PRIME
            assert int(gradient @ matching) % PRIME == 0


def local_generator_evaluations(directions):
    determinants = triple_determinants(test_vectors())
    gradients = {six: gradient_values(six, determinants) for six in SIX_SETS}
    columns = []
    for six, six_gradient in gradients.items():
        complement = tuple(label for label in LABELS if label not in six)
        matching = matching_values(six, directions)
        local_relation = six_gradient @ matching % PRIME
        columns.append(local_relation[:, None] * gradients[complement] % PRIME)
    return np.concatenate(columns, axis=1)


def joint_harmonic_nullity(labels, degree, directions):
    """Nullity of the two squarefree lowering operators at these rows."""
    domain = tuple(combinations(labels, degree))
    codomain = tuple(combinations(labels, degree - 1))
    row_index = {subset: j for j, subset in enumerate(codomain)}
    matrix = np.zeros((2 * len(codomain), len(domain)), dtype=np.int64)
    for column, subset in enumerate(domain):
        for label in subset:
            row = row_index[tuple(j for j in subset if j != label)]
            matrix[row, column] = 1
            matrix[len(codomain) + row, column] = directions[label]
    return len(domain) - modular_rank(matrix)


def check_six_row_core_valuation():
    # The Astra four-row cofactor is F*M with these F edges and row cuts.
    f_edges = ((0, 1), (4, 5), (2, 5), (0, 5), (3, 4))
    row_cuts = ((0, 1), (3, 4), (1, 2), (2, 3))
    representatives = {
        0: (0, 1, 2, 3, 4, 5),
        1: (0, 1, 2, 3, 4, 5),
        2: (0, 2, 1, 3, 4, 5),
        3: (0, 1, 2, 3, 4, 5),
        4: (0, 1, 2, 4, 3, 5),
        5: (0, 1, 2, 3, 4, 5),
        6: (0, 1, 2, 3, 4, 5),
    }
    expected = ((0, 0, 0), (0, 0, 0), (0, 0, 0), (1, 1, 0),
                (1, 2, 1), (2, 4, 2), (5, 8, 3))
    for size, permutation in representatives.items():
        core = set(range(size))
        f = sum(permutation[i] in core and permutation[j] in core
                for i, j in f_edges)
        d = sum(max(0, len(core - {permutation[i], permutation[j]}) - 2)
                for i, j in row_cuts)
        g = max(0, size - 3)
        assert (f, d, g) == expected[size]
        assert f - d + g == 0


def twelve_row_pure_core_initial_matrix(size):
    """Leading normalized coherent matrix for T={0,...,size-1}."""
    core = set(range(size))
    # In projective local coordinates, core rows have slope t*a_i and
    # noncore rows have slope b_i.  The chosen residues are nonzero and
    # distinct where required modulo 101.
    a = [i + 1 for i in LABELS]
    b = [i + 17 for i in LABELS]
    order, leading = {}, {}
    for i, j in combinations(LABELS, 2):
        if i in core and j in core:
            order[(i, j)], leading[(i, j)] = 1, a[j] - a[i]
        elif i in core:
            order[(i, j)], leading[(i, j)] = 0, b[j]
        elif j in core:
            order[(i, j)], leading[(i, j)] = 0, -b[i]
        else:
            order[(i, j)], leading[(i, j)] = 0, b[j] - b[i]

    result = np.zeros((len(CUTS), len(TRIANGLE_COLUMNS)), dtype=np.int64)
    for row, cut in enumerate(CUTS):
        inside = set(cut)
        mandatory = max(0, len(core - inside) - 4)
        for column, triangles in enumerate(TRIANGLE_COLUMNS):
            exponent, coefficient = 0, 1
            for triangle in triangles:
                hits = [j for j, label in enumerate(triangle) if label in inside]
                if len(hits) != 1:
                    coefficient = 0
                    break
                pair = tuple(label for label in triangle if label not in inside)
                exponent += order[pair]
                orientation = -1 if hits[0] == 1 else 1
                coefficient = coefficient * orientation * leading[pair] % PRIME
            if exponent == mandatory:
                result[row, column] = coefficient
    return result


def check_twelve_row_local_core_span():
    """Optional ~25-second initial-syzygy rank audit on all core sizes."""
    determinants = triple_determinants(test_vectors())
    gradients = {six: gradient_values(six, determinants) for six in SIX_SETS}
    a = [i + 1 for i in LABELS]
    b = [i + 17 for i in LABELS]
    for size in range(13):
        core = set(range(size))
        columns = []
        for six, six_gradient in gradients.items():
            mandatory = max(0, len(core.intersection(six)) - 3)
            matching_initial = []
            for matching in MATCHINGS:
                exponent, leading = 0, 1
                for i, j in matching:
                    left, right = six[i], six[j]
                    if left in core and right in core:
                        exponent += 1
                        bracket = a[right] - a[left]
                    elif left in core:
                        bracket = b[right]
                    elif right in core:
                        bracket = -b[left]
                    else:
                        bracket = b[right] - b[left]
                    leading = leading * bracket % PRIME
                matching_initial.append(leading if exponent == mandatory else 0)
            local_relation = six_gradient @ np.array(matching_initial) % PRIME
            complement = tuple(label for label in LABELS if label not in six)
            columns.append(local_relation[:, None] * gradients[complement] % PRIME)
        evaluation_matrix = np.concatenate(columns, axis=1)
        assert modular_rank(evaluation_matrix.T) == 341


def main():
    assert len(TABLEAUX) == 462 and len(CUTS) == 495
    directions = list(LABELS)
    check_local_six_row_incidence(directions)
    determinants = triple_determinants(test_vectors())
    tableau_evaluations = np.column_stack([
        np.prod([determinants[triangle] for triangle in columns], axis=0) % PRIME
        for columns in TRIANGLE_COLUMNS
    ])
    assert tableau_evaluations.shape == (462, 462)
    assert modular_rank(tableau_evaluations) == 462
    coherent = coherent_matrix(directions)
    local = local_generator_evaluations(directions)
    assert local.shape == (462, 4620)
    coherent_rank = modular_rank(coherent)
    local_rank = modular_rank(local.T)
    assert (coherent_rank, local_rank) == (121, 341)
    assert coherent_rank + local_rank == 462
    assert joint_harmonic_nullity(LABELS, 4, directions) == 121
    assert joint_harmonic_nullity(LABELS[1:], 4, directions) == 55
    assert joint_harmonic_nullity(LABELS[1:], 3, directions) == 66

    changed = list(directions)
    changed[0] = 13
    at_13 = coherent_matrix(changed)
    slope = (at_13 - coherent) * pow(13, -1, PRIME) % PRIME
    changed[0] = 29
    assert np.array_equal(coherent_matrix(changed),
                          (coherent + 29 * slope) % PRIME)
    assert (modular_rank(slope), modular_rank(np.concatenate((coherent, slope)))) == (121, 176)
    check_six_row_core_valuation()
    assert [modular_rank(twelve_row_pure_core_initial_matrix(size))
            for size in range(13)] == [121] * 13
    if "--extended-core" in sys.argv[1:]:
        check_twelve_row_local_core_span()
    print("m=12: 462 tableau vectors independent; coherent rank 121,")
    print("      localized-six-row span rank 341, sum 462 mod 101.")
    print("Joint harmonic nullities: (12,4)=121, (11,4)=55, (11,3)=66 mod 101.")
    print("One-row pencil: constant rank 121, slope rank 121, stacked rank 176 mod 101.")
    print("m=6: pure-core normalized maximal-minor valuation zero for all core sizes 0..6.")
    print("m=12: pure-core normalized initial rank 121 for every core size 0..12 at explicit residues.")
    if "--extended-core" in sys.argv[1:]:
        print("m=12: normalized local-six-row initial span rank 341 for every core size 0..12.")


if __name__ == "__main__":
    main()
