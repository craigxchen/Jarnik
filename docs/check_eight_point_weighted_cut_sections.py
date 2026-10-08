#!/usr/bin/env python3
"""Exact sparse checks for weighted clean-cut orders in degree two.

The only finite-field step is a complete intersection-lattice enumeration
inside the canonical 2- and 5-dimensional row-weight spaces.  Reduction can
only delete characteristic-zero support, so its lower bound applies over Q.
"""

from collections import defaultdict, deque
from itertools import combinations, permutations


PRIME = 1009
LABELS = tuple(range(8))
CUTS = [
    cut
    for size in range(1, 5)
    for cut in combinations(LABELS, size)
    if size < 4 or 0 in cut
]


def permutation_sign(values):
    inversions = sum(
        values[i] > values[j]
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def phase_plucker(triple):
    """R_I=(prod_{i in I} z_i^2) Vandermonde(z_j:j not in I)."""
    complement = [i for i in LABELS if i not in triple]
    answer = {}
    for assignment in permutations(range(5)):
        exponent = [0] * 8
        for i in triple:
            exponent[i] = 2
        for i, power in zip(complement, assignment):
            exponent[i] = power
        answer[tuple(exponent)] = permutation_sign(assignment)
    return answer


def multiply(left, right):
    answer = defaultdict(int)
    for a, coefficient_a in left.items():
        for b, coefficient_b in right.items():
            answer[tuple(x + y for x, y in zip(a, b))] += coefficient_a * coefficient_b
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def support_width(polynomial):
    maxima = [-10**9] * len(CUTS)
    for exponent in polynomial:
        centered = tuple(power - 4 for power in exponent)
        for position, cut in enumerate(CUTS):
            maxima[position] = max(
                maxima[position], sum(centered[i] for i in cut)
            )
    return sum(maxima)


def rref(rows, columns):
    rows = [list(value % PRIME for value in row) for row in rows if any(row)]
    answer = []
    for column in range(columns):
        pivot = next((i for i, row in enumerate(rows) if row[column]), None)
        if pivot is None:
            continue
        row = rows.pop(pivot)
        inverse = pow(row[column], PRIME - 2, PRIME)
        row = [value * inverse % PRIME for value in row]
        for collection in (rows, answer):
            for i, other in enumerate(collection):
                if other[column]:
                    multiple = other[column]
                    collection[i] = [
                        (x - multiple * y) % PRIME for x, y in zip(other, row)
                    ]
        answer.append(row)
    return tuple(tuple(row) for row in answer)


def nullspace(rows, columns):
    reduced = rref(rows, columns)
    pivots = [next(i for i, value in enumerate(row) if value) for row in reduced]
    free = [i for i in range(columns) if i not in pivots]
    answer = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in zip(reduced, pivots):
            vector[pivot] = -row[free_column] % PRIME
        answer.append(tuple(vector))
    return tuple(answer)


def canonical_subspace(vectors, ambient_dimension):
    return rref(vectors, ambient_dimension)


def restricted_row(row, subspace):
    return tuple(
        sum(row[i] * basis_vector[i] for i in range(len(row))) % PRIME
        for basis_vector in subspace
    )


def intersect_kernel(subspace, rows, ambient_dimension):
    equations = [restricted_row(row, subspace) for row in rows]
    kernel = nullspace(equations, len(subspace))
    vectors = [
        tuple(
            sum(subspace[j][i] * parameter[j] for j in range(len(subspace)))
            % PRIME
            for i in range(ambient_dimension)
        )
        for parameter in kernel
    ]
    return canonical_subspace(vectors, ambient_dimension)


def canonical_line(line):
    assert len(line) == 1
    vector = line[0]
    pivot = next(value for value in vector if value)
    inverse = pow(pivot, PRIME - 2, PRIME)
    return tuple(value * inverse % PRIME for value in vector)


def filtration_lattice(polynomials):
    """Enumerate every locus on which at least one cut support drops."""
    dimension = len(polynomials)
    exponents = set().union(*(polynomial.keys() for polynomial in polynomials))
    coefficients = {
        exponent: tuple(
            polynomial.get(exponent, 0) % PRIME for polynomial in polynomials
        )
        for exponent in exponents
    }
    assert len(rref(coefficients.values(), dimension)) == dimension

    cut_levels = []
    for cut in CUTS:
        levels = defaultdict(set)
        for exponent, row in coefficients.items():
            weight = sum(exponent[i] - 4 for i in cut)
            levels[weight].add(row)
        cut_levels.append(
            [
                (weight, rref(rows, dimension))
                for weight, rows in sorted(levels.items(), reverse=True)
            ]
        )

    whole_space = canonical_subspace(
        [tuple(int(i == j) for i in range(dimension)) for j in range(dimension)],
        dimension,
    )
    queue = deque([whole_space])
    seen = {whole_space}
    widths = {}
    while queue:
        subspace = queue.popleft()
        width = 0
        children = []
        for levels in cut_levels:
            for weight, rows in levels:
                if any(any(restricted_row(row, subspace)) for row in rows):
                    width += weight
                    child = intersect_kernel(subspace, rows, dimension)
                    if child and len(child) < len(subspace):
                        children.append(child)
                    break
        widths[subspace] = width
        for child in children:
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return widths


def raw_cut_order(triple, cut):
    intersection = len(set(triple) & set(cut))
    remaining = len(cut) - intersection
    return 2 * intersection + remaining * (remaining - 1) // 2


def normalized_cut_order(pair, cut):
    minimum = (0, 1, 3, 5)[len(cut) - 1]
    return sum(raw_cut_order(triple, cut) for triple in pair) - 2 * minimum


def check_all_coordinate_products():
    triples = list(combinations(LABELS, 3))
    for position, left in enumerate(triples):
        for right in triples[position:]:
            total = sum(normalized_cut_order((left, right), cut) for cut in CUTS)
            assert total == 106


def check_reciprocity():
    for triple in combinations(LABELS, 3):
        polynomial = phase_plucker(triple)
        assert all(
            polynomial.get(tuple(4 - power for power in exponent)) == coefficient
            for exponent, coefficient in polynomial.items()
        )


def weight_space_polynomials(specifications):
    needed = {triple for pair in specifications for triple in pair}
    pluckers = {triple: phase_plucker(triple) for triple in needed}
    return [multiply(pluckers[left], pluckers[right]) for left, right in specifications]


def check_weight_spaces():
    repeated_label_basis = [
        ((0, 1, 2), (0, 3, 4)),
        ((0, 1, 3), (0, 2, 4)),
    ]
    repeated = weight_space_polynomials(repeated_label_basis)
    repeated_lattice = filtration_lattice(repeated)
    assert len(repeated_lattice) == 4
    assert min(repeated_lattice.values()) == 640
    assert sorted(len(space) for space in repeated_lattice) == [1, 1, 1, 2]
    repeated_equality = {
        canonical_line(space)
        for space, width in repeated_lattice.items()
        if len(space) == 1 and width == 640
    }
    assert repeated_equality == {
        (1, 0),
        (0, 1),
        (1, PRIME - 1),
    }

    disjoint_basis = [
        ((0, 1, 2), (3, 4, 5)),
        ((0, 1, 3), (2, 4, 5)),
        ((0, 1, 4), (2, 3, 5)),
        ((0, 2, 3), (1, 4, 5)),
        ((0, 2, 4), (1, 3, 5)),
    ]
    disjoint = weight_space_polynomials(disjoint_basis)
    assert all(support_width(polynomial) == 640 for polynomial in disjoint)
    assert max(abs(value) for polynomial in disjoint for value in polynomial.values()) == 2
    disjoint_lattice = filtration_lattice(disjoint)
    histogram = defaultdict(int)
    for space in disjoint_lattice:
        histogram[len(space)] += 1
    assert dict(histogram) == {5: 1, 3: 15, 1: 25, 2: 45}
    assert min(disjoint_lattice.values()) == 640
    disjoint_equality = {
        canonical_line(space)
        for space, width in disjoint_lattice.items()
        if len(space) == 1 and width == 640
    }
    expected_monomial_lines = {
        (1, 0, 0, 0, 0),
        (0, 1, 0, 0, 0),
        (0, 0, 1, 0, 0),
        (0, 0, 0, 1, 0),
        (0, 0, 0, 0, 1),
        (1, PRIME - 1, 1, 0, 0),
        (1, 0, 0, 1, PRIME - 1),
        (1, PRIME - 1, 1, 1, PRIME - 1),
        (1, 0, 1, 0, PRIME - 1),
        (1, PRIME - 1, 0, 1, 0),
    }
    assert disjoint_equality == expected_monomial_lines


def main():
    assert len(CUTS) == 127
    check_reciprocity()
    check_all_coordinate_products()
    check_weight_spaces()
    print(
        "PASS: exact cut filtrations give width >= 640 in every canonical "
        "degree-two row-weight space."
    )


if __name__ == "__main__":
    main()
