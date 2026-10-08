"""Exact matching lifts and integral numerical-restriction index checks."""

from itertools import combinations
from math import gcd, prod

from check_all_cut_invariant_relation_hierarchy import (
    add, basis, harmonic_basis, restrict_graph,
)


def graph_value(graph, directions):
    return prod(directions[i][0]*directions[j][1]
                - directions[i][1]*directions[j][0] for i, j in graph)


def lift_graph(m, k, pairs):
    inside_count = m//2-k
    inside = list(range(inside_count))
    outside = list(range(inside_count, m))
    graph = []
    used = set()
    for s, (a, b) in zip(inside, pairs):
        a, b = outside[a], outside[b]
        graph.extend(((s, a), (a, b), (b, s)))
        used.update((a, b))
    unused = [a for a in outside if a not in used]
    assert len(unused) == len(inside[2*k:]) == m//2-3*k
    for a, b in zip(inside[2*k:], unused):
        graph.extend(((a, b), (a, b)))
    return graph, set(inside), outside


def check_lifts():
    count = 0
    for m in (6, 8, 10, 12):
        directions = [(j+1, 1) for j in range(m)]
        for k in range(1, m//6+1):
            for _, pairs, polynomial in harmonic_basis(m//2+k, 2*k):
                graph, inside, outside = lift_graph(m, k, pairs)
                assert all(sum(j in edge for edge in graph) == 2 for j in range(m))
                expected = {}
                for mask, coefficient in polynomial.items():
                    exponents = [0]*m
                    for j, label in enumerate(outside):
                        exponents[label] = int(bool(mask & (1 << j)))
                    expected[tuple(exponents)] = coefficient
                assert restrict_graph(graph, inside, m) == expected
                # Every larger inside set contains an internal component edge.
                for selected in combinations(range(m), len(inside)+1):
                    chosen = set(selected)
                    assert any(a in chosen and b in chosen for a, b in graph)
                outside_first = outside[pairs[0][0]]
                swap = {1: outside_first, outside_first: 1}
                witness = [(swap.get(a, a), swap.get(b, b)) for a, b in graph]
                assert not restrict_graph(witness, inside, m)
                assert all(sum(j in edge for edge in witness) == 2 for j in range(m))
                for selected in combinations(range(m), len(inside)+1):
                    chosen = set(selected)
                    assert any(a in chosen and b in chosen for a, b in witness)
                fQ = graph_value(graph, directions)
                fZ = graph_value(witness, directions)
                assert fQ != 0 and fZ != 0
                count += 1
    return count


def extended_gcd(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        quotient = old_r//r
        old_r, r = r, old_r-quotient*r
        old_s, s = s, old_s-quotient*s
        old_t, t = t, old_t-quotient*t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def kernel_basis(f):
    size = len(f)
    columns = [[int(i == j) for i in range(size)] for j in range(size)]
    values = list(f)
    for j in range(1, size):
        a, b = values[0], values[j]
        if a == b == 0:
            continue
        g, s, t = extended_gcd(a, b)
        left, right = columns[0], columns[j]
        columns[0] = [s*x+t*y for x, y in zip(left, right)]
        columns[j] = [(-b//g)*x+(a//g)*y for x, y in zip(left, right)]
        values[0], values[j] = g, 0
    assert values[0] == gcd(*f)
    assert all(sum(a*b for a, b in zip(f, column)) == 0 for column in columns[1:])
    assert abs(determinant([[columns[j][i] for j in range(size)]
                            for i in range(size)])) == 1
    return columns[1:]


def determinant(matrix):
    if not matrix:
        return 1
    if len(matrix) == 1:
        return matrix[0][0]
    return sum((-1)**j*entry*determinant([row[:j]+row[j+1:] for row in matrix[1:]])
               for j, entry in enumerate(matrix[0]) if entry)


def minors_gcd(columns, rows, order):
    answer = 0
    for selected_rows in combinations(range(rows), order):
        for selected_columns in combinations(range(len(columns)), order):
            value = determinant([[columns[j][i] for j in selected_columns]
                                 for i in selected_rows])
            answer = gcd(answer, value)
    return answer


def check_indices():
    count = 0
    for size in range(3, 7):
        for target_rank in range(1, size):
            for seed in range(1, 25):
                f = [((seed+3)*(j+2)+j*j) % 17-8 for j in range(size)]
                if not any(f[target_rank:]):
                    f[-1] = 6
                if seed % 3 == 0:
                    f = [6*v for v in f]
                columns = kernel_basis(f)
                projected = [column[:target_rank] for column in columns]
                a, b = gcd(*f[target_rank:]), gcd(*f)
                assert a > 0 and b > 0 and a % b == 0
                assert minors_gcd(projected, target_rank, target_rank) == a//b
                # All smaller determinantal divisors are 1: the quotient is cyclic.
                for order in range(1, target_rank):
                    assert minors_gcd(projected, target_rank, order) == 1
                # Also verify the index relative to a nonprimitive ambient image I.
                scales = [j+2 for j in range(target_rank)]
                scaled = [[value*factor for value, factor in zip(column, scales)]
                          for column in projected]
                assert minors_gcd(scaled, target_rank, target_rank) == a//b*prod(scales)
                count += 1
    return count


def split_columns(matrix):
    """Unimodular column operations, with exact Bezout certificates."""
    reduced = [row[:] for row in matrix]
    row_count, column_count = len(matrix), len(matrix[0])
    transform = [[int(i == j) for j in range(column_count)]
                 for i in range(column_count)]
    for i in range(row_count):
        pivot = next(j for j in range(i, column_count) if reduced[i][j])
        for target in (reduced, transform):
            for row in target:
                row[i], row[pivot] = row[pivot], row[i]
        for j in range(i+1, column_count):
            a, b = reduced[i][i], reduced[i][j]
            if not b:
                continue
            g, s, t = extended_gcd(a, b)
            assert s*a+t*b == g and g > 0
            # The following two-column transformation has determinant one.
            for target in (reduced, transform):
                for row in target:
                    left, right = row[i], row[j]
                    row[i] = s*left+t*right
                    row[j] = (-b//g)*left+(a//g)*right
    assert all(reduced[i][j] == sum(matrix[i][a]*transform[a][j]
                                  for a in range(column_count))
               for i in range(row_count) for j in range(column_count))
    assert all(not reduced[i][j] for i in range(row_count)
               for j in range(row_count, column_count))
    return reduced, transform


def polynomial_multiply(a, b):
    result = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def check_realized_spanning_family():
    graphs = basis(6)
    assert len(graphs) == 15
    balanced = []
    for selected in combinations(range(6), 3):
        if 0 not in selected:
            continue
        inside = set(selected)
        balanced.append([prod(0 if (a in inside) == (b in inside)
                              else (1 if a in inside else -1) for a, b in graph)
                         for graph in graphs])
    _, transform = split_columns(balanced)
    kernel_columns = [[row[j] for row in transform] for j in range(10, 15)]
    restrictions = [restrict_graph(graph, {0, 1}, 6) for graph in graphs]
    raw_restriction = []
    for a, b in ((2, 4), (2, 3)):
        monomial = tuple(int(j in (a, b)) for j in range(6))
        raw_restriction.append([p.get(monomial, 0) for p in restrictions])
    restriction = [[sum(a*b for a, b in zip(row, column))
                    for column in kernel_columns] for row in raw_restriction]
    reduced, splitting = split_columns(restriction)
    assert reduced == [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0]]
    columns = [[sum(kernel_columns[j][a]*splitting[j][b] for j in range(5))
                for a in range(15)] for b in range(5)]
    terms = [{5: 1, 10: -1, 12: 1}, {5: 1, 8: -1, 10: -1, 13: 1},
             {3: 1}, {7: 1}, {5: 1, 6: -1, 8: 1}]
    assert columns == [[term.get(j, 0) for j in range(15)] for term in terms]
    expected_restrictions = [
        restrict_graph([(3, 2), (5, 4)], set(), 6),
        restrict_graph([(4, 2), (5, 3)], set(), 6), {}, {}, {},
    ]
    for column, expected in zip(columns, expected_restrictions):
        polynomial = {}
        for coefficient, restricted in zip(column, restrictions):
            polynomial = add(polynomial, restricted, coefficient)
        assert polynomial == expected
    # Exact graph polynomials for directions (s,1),(0,1),...,(4,1).
    graph_polynomials = []
    for graph in graphs:
        polynomial = [1]
        for a, b in graph:
            left = [0, 1] if a == 0 else [a-1, 0]
            right = [0, 1] if b == 0 else [b-1, 0]
            factor = [left[j]-right[j] for j in range(2)]
            polynomial = polynomial_multiply(polynomial, factor)
        assert not any(polynomial[3:])
        graph_polynomials.append(polynomial[:3])
    polynomials = [[sum(column[j]*graph_polynomials[j][degree] for j in range(15))
                    for degree in range(3)] for column in columns]
    assert polynomials == [[24, -38, 14], [96, -74, 14], [0, -2, 2],
                           [0, -24, 12], [0, -56, 20]]
    u, constant_kernel = [-4, 1, 3, 3, 0], [0, 0, 8, -3, 1]
    for relation in (u, constant_kernel):
        assert all(sum(c*p[d] for c, p in zip(relation, polynomials)) == 0
                   for d in range(3))
    # v_s=(s,0,12-7s,0,0) is an identity of polynomial evaluations.
    left = polynomial_multiply([0, 1], polynomials[0])
    right = polynomial_multiply([12, -7], polynomials[2])
    assert all(a+b == 0 for a, b in zip(left, right))
    for s in range(5, 125, 6):
        f = [a+b*s+c*s*s for a, b, c in polynomials]
        assert gcd(*f) == 4 and gcd(*f[2:]) == 4*s
        assert f[0]//4 % s == 6 % s and f[1]//4 % s == 24 % s
        numerical_kernel = kernel_basis(f)
        projected = [column[:2] for column in numerical_kernel]
        assert minors_gcd(projected, 2, 2) == s
        assert all((x+4*y) % s == 0 for x, y in projected)
        assert determinant([[u[0], s], [u[1], 0]]) == -s
        assert all((x+4*y == 0) for x in range(-(s-1)//5, (s-1)//5+1)
                   for y in range(-(s-1)//5, (s-1)//5+1)
                   if 5*abs(x) < s and 5*abs(y) < s and (x+4*y) % s == 0)
    return 20


if __name__ == '__main__':
    print('Passed:', check_lifts(), 'integral matching lifts and actual kernel witnesses;',
          check_indices(), 'full integer relation-image index and cyclic quotient cases;',
          check_realized_spanning_family(), 'realized spanning-height specializations '
          'with an exact saturated-basis and polynomial certificate.')
