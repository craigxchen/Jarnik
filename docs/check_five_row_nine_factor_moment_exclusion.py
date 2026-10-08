"""Zero-edge certificates and repeated-threshold reduction for nine factors."""

from itertools import combinations

from check_k5_real_moment_exclusion import (
    CERTIFICATES, EDGES, PATTERNS, RAYS, UNITS, WORDS, ZERO,
    add, admissible_words, canonical, linear, multiply, rank,
)


def main():
    critical = (1, 1, -1, -1, 1)
    for prefix in (1, -1):
        for sign in (1, -1):
            assert (prefix,) + tuple(sign * x for x in critical) in PATTERNS
    words = admissible_words()
    assert len(words) == 48
    assert sorted({canonical(word) for word in words}) == list(WORDS)

    # The length-nine combinatorial boundary exists; parity completes it
    # uniquely to all ten K5 columns. Its moment magnitudes are excluded below.
    rows9 = [tuple(int(i in edge) for edge in EDGES[1:]) for i in range(5)]
    rows10 = [row + (sum(row) % 2,) for row in rows9]
    assert min(sum(a != b for a, b in zip(rows9[i], rows9[j]))
               for i, j in EDGES) == 5
    assert all(sum(a != b for a, b in zip(rows10[i], rows10[j])) == 6
               for i, j in EDGES)
    columns = [tuple(i for i in range(5) if rows10[i][j]) for j in range(10)]
    assert sorted(columns) == list(EDGES)
    incidence = [[int((i in edge) != (j in edge)) for edge in EDGES]
                 for i, j in EDGES]
    assert rank(incidence) == 10
    assert all(sum(row) == 6 for row in incidence)

    # A repeated positive coefficient supplies nested threshold supports
    # before minority orientation. Enumerate every possible nonconstant
    # support surviving the parity equality and every strict nested pair.
    vertices = frozenset(range(5))
    supports = [frozenset(s) for size in (2, 3)
                for s in combinations(range(5), size)]
    minority = lambda s: s if len(s) == 2 else vertices - s
    nested = [(a, b) for a in supports for b in supports if a < b]
    assert len(nested) == 30
    assert all(len(a) == 2 and len(b) == 3 and
               minority(a).isdisjoint(minority(b)) for a, b in nested)
    assert not any(a < b < c for a in supports for b in supports for c in supports)

    expected = ((0, 31, 3), (1, 48, 108), (2, 33, 3), (2, 45, 3870))
    for number, (word, rays, certificate, target) in enumerate(
            zip(WORDS, RAYS, CERTIFICATES, expected)):
        coefficients = [[sign * (int(vertex in EDGES[e]) - int(0 in EDGES[e]))
                         for e, sign in word] for vertex in range(1, 5)]
        kernel_matrix = [[sum(row[j:]) for j in range(10)] for row in coefficients]
        assert rank(kernel_matrix) == 4 and rank(rays) == 6
        assert all(sum(a * b for a, b in zip(row, ray)) == 0
                   for row in kernel_matrix for ray in rays)
        inverse = [next(i for i in range(10)
                        if tuple(ray[i] for ray in rays) == UNITS[j])
                   for j in range(6)]
        zero_parameter = inverse.index(0)
        assert all(i > 0 for j, i in enumerate(inverse) if j != zero_parameter)
        gaps = [linear([ray[i] for ray in rays]) for i in range(10)]
        magnitudes, partial = [], {}
        for gap in gaps:
            partial = add(partial, gap)
            magnitudes.append(partial)
        cubes = [multiply(multiply(x, x), x) for x in magnitudes]
        equations = []
        for row in coefficients:
            f = {}
            for weight, cubic in zip(row, cubes):
                f = add(f, cubic, weight)
            equations.append(f)
        polynomial = {}
        for weights, f in zip(certificate, equations):
            multiplier = {ZERO: weights[0]} if len(weights) == 1 else linear(weights)
            polynomial = add(polynomial, multiply(multiplier, f))
        boundary = {monomial: value for monomial, value in polynomial.items()
                    if monomial[zero_parameter] == 0}
        assert boundary and all(value > 0 for value in boundary.values())
        assert (zero_parameter, len(boundary), min(boundary.values())) == target
        print(f'PASS orbit {number}: zero lambda {zero_parameter}, '
              f'{len(boundary)} positive boundary terms, minimum {min(boundary.values())}')
    print('PASS: complete K5 zero-edge exclusion, length-nine parity reduction, '
          'and all 30 repeated-threshold nested pairs')


if __name__ == '__main__':
    main()
