"""Exact polynomial checks for the all-cut hierarchy; standard library only.

The arithmetic prime-power estimates are proved in the companion note.
These checks certify the matching bases, integral conversion, commutator,
finite low-row restriction ranks, and higher-depth triangle witnesses.
"""

from itertools import combinations, product
from math import comb, prod


def basis(m):
    answer = []
    for top_counts in product(range(3), repeat=m):
        if sum(top_counts) != m:
            continue
        top = [i for i, count in enumerate(top_counts) for _ in range(count)]
        bottom = [i for i, count in enumerate(top_counts) for _ in range(2 - count)]
        if all(a < b for a, b in zip(top, bottom)):
            answer.append(list(zip(top, bottom)))
    return answer


def dimension(degrees):
    coefficients = [1]
    for degree in degrees:
        updated = [0] * (len(coefficients) + degree)
        for i, coefficient in enumerate(coefficients):
            for j in range(degree + 1):
                updated[i + j] += coefficient
        coefficients = updated
    midpoint = sum(degrees) // 2
    return coefficients[midpoint] - coefficients[midpoint - 1]


def add(left, right, scale=1):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, 0) + scale * coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def down(poly, n):
    out = {}
    for mask, coefficient in poly.items():
        for j in range(n):
            if mask & (1 << j):
                out = add(out, {mask ^ (1 << j): coefficient})
    return out


def up(poly, n):
    out = {}
    for mask, coefficient in poly.items():
        for j in range(n):
            if not mask & (1 << j):
                out = add(out, {mask | (1 << j): coefficient})
    return out


def matching_polynomial(pairs):
    out = {0: 1}
    for a, b in pairs:
        updated = {}
        for mask, coefficient in out.items():
            assert not mask & ((1 << a) | (1 << b))
            updated = add(updated, {mask | (1 << b): coefficient})
            updated = add(updated, {mask | (1 << a): -coefficient})
        out = updated
    return out


def harmonic_basis(n, d):
    answer = []
    for bottom in combinations(range(n), d):
        if any(b < 2 * j + 1 for j, b in enumerate(bottom)):
            continue
        unused = set(range(n)) - set(bottom)
        pairs = []
        for b in bottom:
            a = min(j for j in unused if j < b)
            unused.remove(a)
            pairs.append((a, b))
        mask = sum(1 << b for b in bottom)
        answer.append((mask, pairs, matching_polynomial(pairs)))
    return sorted(answer)


def check_bases():
    for n in range(1, 13):
        for d in range(n + 1):
            for indices in combinations(range(n), d):
                mask = sum(1 << j for j in indices)
                commutator = add(down(up({mask: 1}, n), n),
                                 up(down({mask: 1}, n), n), -1)
                target = {mask: n - 2 * d} if n != 2 * d else {}
                assert commutator == target
            if d > n // 2:
                continue
            polynomials = harmonic_basis(n, d)
            expected = comb(n, d) - (comb(n, d - 1) if d else 0)
            assert len(polynomials) == expected
            for mask, pairs, polynomial in polynomials:
                assert len(pairs) == d and len({j for pair in pairs for j in pair}) == 2 * d
                assert polynomial[mask] == 1 and max(polynomial) == mask
                assert all(abs(c) == 1 for c in polynomial.values())
                assert not down(polynomial, n)
            # Recover arbitrary integer coordinates by triangular elimination.
            original = [(j % 7) - 3 for j in range(expected)]
            value = {}
            for coefficient, (_, _, polynomial) in zip(original, polynomials):
                value = add(value, polynomial, coefficient)
            recovered = [0] * expected
            residue = dict(value)
            for j in reversed(range(expected)):
                mask, _, polynomial = polynomials[j]
                recovered[j] = residue.get(mask, 0)
                residue = add(residue, polynomial, -recovered[j])
            assert not residue and recovered == original
            assert sum(abs(c) for c in recovered) <= 2 ** (expected - 1) * sum(abs(c) for c in value.values())
    print("Squarefree D/U identity, ballot matching bases, and integral conversion pass for n<=12.")


def restrict_graph(graph, inside, m):
    # Inside=(0,1), outside=(1,z). Bracket convention x_i*y_j-y_i*x_j.
    polynomial = {(0,) * m: 1}
    for i, j in graph:
        if i in inside and j in inside:
            return {}
        if i in inside:
            polynomial = {key: -c for key, c in polynomial.items()}
        elif j not in inside:
            updated = {}
            for exponents, coefficient in polynomial.items():
                for label, sign in ((j, 1), (i, -1)):
                    new = list(exponents)
                    new[label] += 1
                    key = tuple(new)
                    updated[key] = updated.get(key, 0) + sign * coefficient
            polynomial = {key: c for key, c in updated.items() if c}
    return polynomial


def modular_rank(rows, columns, prime=1000003):
    pivots = {}
    for row in rows:
        row = [value % prime for value in row]
        for j in sorted(pivots):
            if row[j]:
                multiplier = row[j]
                row = [(a - multiplier * b) % prime for a, b in zip(row, pivots[j])]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            inverse = pow(row[pivot], prime - 2, prime)
            pivots[pivot] = [(value * inverse) % prime for value in row]
            if len(pivots) == columns:
                break
    return len(pivots)


def restriction_rows(graphs, m, size):
    for selected in combinations(range(m), size):
        restrictions = [restrict_graph(graph, set(selected), m) for graph in graphs]
        monomials = set().union(*(set(p) for p in restrictions))
        for monomial in sorted(monomials):
            yield [p.get(monomial, 0) for p in restrictions]


def check_low_row_ranks():
    for m in (6, 8):
        graphs = basis(m)
        assert len(graphs) == dimension([2] * m)
        balanced_rank = modular_rank(restriction_rows(graphs, m, m // 2), len(graphs))
        full_rank = modular_rank(restriction_rows(graphs, m, m // 2 - 1), len(graphs))
        assert balanced_rank == comb(m, m // 2) // 2
        assert full_rank == len(graphs)
        print(f"m={m}: balanced rank {balanced_rank}; first-smaller joint rank {full_rank} (full).")


def check_triangle_witnesses():
    for m in range(6, 21, 2):
        q = m // 2
        for k in range(1, m // 6 + 1):
            graph = []
            selected = set()
            for start in range(0, 6 * k, 3):
                graph += [(start, start + 1), (start + 1, start + 2), (start + 2, start)]
                selected.add(start)
            for start in range(6 * k, m, 2):
                graph += [(start, start + 1)] * 2
                selected.add(start)
            assert len(selected) == q - k
            polynomial = restrict_graph(graph, selected, m)
            assert polynomial
            assert all(sum(exponents) == 2 * k and max(exponents) == 1 for exponents in polynomial)
            # An independent set has <=1 vertex in each triangle or doubled edge.
            assert 2 * k + (m - 6 * k) // 2 == q - k
            for label in set(range(m)) - selected:
                assert not restrict_graph(graph, selected | {label}, m)
    print("Nonzero triangle-product witnesses at every allowable depth pass through m=20.")


def check_denominator_clearing():
    # Rational slopes Y/P are checked after multiplying by all outside P's.
    from fractions import Fraction
    for n, d in ((4, 2), (6, 2), (8, 4), (9, 4), (12, 6)):
        p = [2 * j + 3 for j in range(n)]
        y = [j * j + 1 for j in range(n)]
        z = [Fraction(b, a) for a, b in zip(p, y)]
        for _, pairs, polynomial in harmonic_basis(n, d):
            value = sum(c * prod(z[j] for j in range(n) if mask & (1 << j))
                        for mask, c in polynomial.items())
            used = {j for pair in pairs for j in pair}
            bracket_value = prod(y[b] * p[a] - y[a] * p[b] for a, b in pairs)
            bracket_value *= prod(p[j] for j in range(n) if j not in used)
            assert prod(p) * value == bracket_value
    for m in range(6, 101, 2):
        for k in range(1, m // 6 + 1):
            n = m // 2 + k
            one_pass = 2 * n + 2 * (n - 4 * k) + 4 * k
            two_pass = one_pass + (2 * m - 1) * one_pass + 2 * m
            assert one_pass == 2 * m
            assert two_pass == 2 * m * (2 * m + 1)
    print("Exact matching denominator clearing and both finite-loss identities pass.")


if __name__ == "__main__":
    check_bases()
    check_low_row_ranks()
    check_triangle_witnesses()
    check_denominator_clearing()
    print("All-cut hierarchy checks passed; no rank-to-zero induction or uniform bound asserted.")
