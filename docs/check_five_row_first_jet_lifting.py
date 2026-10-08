"""Exact first-jet lattice and common-Q stress lifting on five rows.

No bound on the coefficient height of the lifted invariant is asserted.
"""

from functools import reduce
from itertools import combinations, combinations_with_replacement
from math import gcd, prod

from check_five_row_gradient_cofactor_stress import (
    MINIMUM, bracket, grad_order, mul as gaussian_mul, norm, power, split_primes,
)
from check_five_row_gradient_discriminant_core_count import (
    determinant, evaluate, independent_graphs, independent_monomials,
    reduce_with_pivots, solve_coefficients,
)
from check_smaller_cut_relation_rank import add, derivative, mul, scale


TRIPLES = list(combinations(range(6), 3))
PRODUCTS = list(combinations_with_replacement(range(6), 2))


def chart(poly):
    out = {}
    for m, c in poly.items():
        if m[0] or m[3]:
            continue
        key = (m[7], m[9])
        out[key] = out.get(key, 0) + c
    return {m: c for m, c in out.items() if c}


def polynomial_determinant(matrix):
    out = {}
    for permutation, sign in (
        ((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
        ((0, 2, 1), -1), ((1, 0, 2), -1), ((2, 1, 0), -1),
    ):
        term = mul(mul(matrix[0][permutation[0]], matrix[1][permutation[1]]),
                   matrix[2][permutation[2]])
        out = add(out, scale(term, sign))
    return out


def polynomial_basis(polys):
    pivots, chosen = {}, []
    for i, poly in enumerate(polys):
        reduced = reduce_with_pivots(poly, pivots)
        if reduced:
            pivot = min(reduced)
            pivots[pivot] = scale(reduced, 1 / reduced[pivot])
            chosen.append((i, poly))
    return chosen


def rank_mod(rows, prime):
    pivots = {}
    for row in rows:
        row = [v % prime for v in row]
        for pivot, previous in sorted(pivots.items()):
            c = row[pivot]
            row = [(x - c * y) % prime for x, y in zip(row, previous)]
        pivot = next((i for i, c in enumerate(row) if c), None)
        if pivot is not None:
            inverse = pow(row[pivot], -1, prime)
            pivots[pivot] = [(c * inverse) % prime for c in row]
    return len(pivots)


def quadratic_product_data():
    basis = independent_graphs(5, 2)
    sections = [chart(poly) for _, poly in basis]
    jet_sections = [polynomial_determinant([
        [sections[c] for c in triple],
        [derivative(sections[c], 0) for c in triple],
        [derivative(sections[c], 1) for c in triple],
    ]) for triple in TRIPLES]
    products = [mul(sections[i], sections[j]) for i, j in PRODUCTS]
    product_basis = polynomial_basis(products)
    jet_basis = polynomial_basis(jet_sections)
    assert len(product_basis) == len(jet_basis) == 16
    assert len(jet_sections) == 20 and all(jet_sections)

    monomials = independent_monomials(product_basis)
    assert determinant([[p.get(m, 0) for _, p in product_basis]
                        for m in monomials]) == -1
    for p in products:
        assert all(c.denominator == 1 for c in solve_coefficients(product_basis, p))
    jet_coefficients = [solve_coefficients(product_basis, p) for p in jet_sections]
    assert all(c.denominator == 1 for row in jet_coefficients for c in row)
    jet_coefficients = [[int(c) for c in row] for row in jet_coefficients]
    selected = [i for i, _ in jet_basis]
    assert abs(determinant([jet_coefficients[i] for i in selected])) == 6
    assert rank_mod(jet_coefficients, 2) == 16
    assert rank_mod(jet_coefficients, 3) == 15
    # The complete jet-section lattice has index exactly3: it divides6,
    # has no factor2, and does have a factor3.
    one, zero = {(0, 0): 1}, {}
    normalized_rows = [(zero, one), (one, zero), (one, one),
                       (one, {(1, 0): 1}), (one, {(0, 1): 1})]
    # These constant frame choices have determinant one on the chart.
    chart_jets = [[scale(chart(derivative(poly, 0)), -1) for _, poly in basis]]
    chart_jets += [[chart(derivative(poly, 2 * i + 1)) for _, poly in basis]
                   for i in range(1, 5)]
    for i, j in combinations(range(5), 2):
        complement = [k for k in range(5) if k not in (i, j)]
        triangle = scale(one, (-1) ** (i + j + 1))
        for a, b in combinations(complement, 2):
            xa, ya = normalized_rows[a]
            xb, yb = normalized_rows[b]
            triangle = mul(triangle, add(mul(xa, yb), scale(mul(ya, xb), -1)))
        for k, triple in enumerate(TRIPLES):
            actual = polynomial_determinant([
                [sections[c] for c in triple],
                [chart_jets[i][c] for c in triple],
                [chart_jets[j][c] for c in triple],
            ])
            assert actual == mul(triangle, jet_sections[k])
    return basis, product_basis, jet_coefficients


def extended_gcd(a, b):
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    g, s, t = extended_gcd(b, a % b)
    return g, t, s - (a // b) * t


def row_reduction(row):
    """Return unimodular columns U with row*U=(gcd(row),0,...)."""
    n = len(row)
    columns = [[int(i == j) for i in range(n)] for j in range(n)]
    values = list(row)
    for j in range(1, n):
        a, b = values[0], values[j]
        if b == 0:
            continue
        g, s, t = extended_gcd(a, b)
        first, other = columns[0], columns[j]
        columns[0] = [s * x + t * y for x, y in zip(first, other)]
        columns[j] = [-(b // g) * x + (a // g) * y for x, y in zip(first, other)]
        values[0], values[j] = g, 0
    if values[0] < 0:
        columns[0] = [-c for c in columns[0]]
        values[0] = -values[0]
    assert [sum(a * b for a, b in zip(row, c)) for c in columns] == values
    return columns, values[0]


def det3(matrix):
    a, b, c = matrix
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def solve_two_rows(rows, target):
    columns, first_gcd = row_reduction(rows[0])
    assert target[0] % first_gcd == 0
    transformed_second = [sum(a * b for a, b in zip(rows[1], c)) for c in columns]
    z0 = target[0] // first_gcd
    remainder = target[1] - transformed_second[0] * z0
    second_columns, second_gcd = row_reduction(transformed_second[1:])
    assert remainder % second_gcd == 0
    z = [z0] + [c * (remainder // second_gcd) for c in second_columns[0]]
    solution = [sum(columns[j][i] * z[j] for j in range(len(z)))
                for i in range(len(z))]
    assert [sum(a * b for a, b in zip(row, solution)) for row in rows] == list(target)
    return solution


def check_source(basis, product_basis, jet_coefficients, correction_exponents):
    primes = split_primes(32)
    blocks = {u: primes[u - 1] for u in range(1, 32)}
    n = {u: norm(blocks[u]) for u in blocks}
    points = []
    for i, exponent in enumerate(correction_exponents):
        p = power(primes[-1], exponent)
        for u in blocks:
            if u & (1 << i):
                p = gaussian_mul(p, blocks[u])
        assert gcd(abs(p[0]), abs(p[1])) == 1
        points.append(p)
    values = sum(points, ())
    evaluations = [evaluate(poly, values) for _, poly in basis]
    gval = reduce(gcd, (abs(v) for v in evaluations))
    all_products = [evaluations[i] * evaluations[j] for i, j in PRODUCTS]
    product_values = [all_products[index] for index, _ in product_basis]
    svalues = [sum(c * v for c, v in zip(row, product_values))
               for row in jet_coefficients]
    s_gcd = reduce(gcd, (abs(v) for v in svalues))
    assert s_gcd % (gval * gval) == 0
    epsilon = s_gcd // (gval * gval)
    assert epsilon in (1, 3)

    jets = []
    for i, (x, y) in enumerate(points):
        gg, u, v = extended_gcd(x, y)
        assert gg == 1
        frame = (-v, u)
        assert bracket(points[i], frame) == 1
        jets.append([frame[0] * evaluate(derivative(poly, 2 * i), values)
                     + frame[1] * evaluate(derivative(poly, 2 * i + 1), values)
                     for _, poly in basis])
    for i, j in combinations(range(5), 2):
        complement = [k for k in range(5) if k not in (i, j)]
        triangle = ((-1) ** (i + j + 1)
                    * prod(bracket(points[a], points[b])
                           for a, b in combinations(complement, 2)))
        for k, triple in enumerate(TRIPLES):
            actual = det3([[row[c] for c in triple]
                           for row in (evaluations, jets[i], jets[j])])
            assert actual == triangle * svalues[k]

    divisors = [prod(n[u] ** grad_order(u, i) for u in blocks) for i in range(5)]
    D0 = prod(n[u] ** max(0, 2 * bin(u).count("1") - 5) for u in blocks)
    L = prod(n[u] ** MINIMUM[bin(u).count("1")] for u in blocks)
    assert prod(divisors) == D0 * L
    assert gval % D0 == 0
    delta = gval // D0
    aminors = {}
    for triple in combinations(range(5), 3):
        aminors[triple] = (prod(divisors[i] for i in triple)
                          * prod(bracket(points[i], points[j])
                                 for i, j in combinations(triple, 2)))
        assert aminors[triple] % L == 0
    agcd = reduce(gcd, (abs(v) for v in aminors.values()))
    rho = agcd // L
    corrections = []
    for i, j in combinations(range(5), 2):
        pair_core = prod(n[u] for u in blocks if u & (1 << i) and u & (1 << j))
        assert abs(bracket(points[i], points[j])) % pair_core == 0
        corrections.append(abs(bracket(points[i], points[j])) // pair_core)
    bproduct = prod(corrections)
    assert bproduct ** 2 % delta == 0
    assert bproduct % rho == 0

    columns, gg = row_reduction(evaluations)
    assert gg == gval
    kernel = columns[1:]
    image = []
    for i in range(5):
        row = [sum(a * b for a, b in zip(jets[i], c)) for c in kernel]
        assert all(v % divisors[i] == 0 for v in row)
        image.append([v // divisors[i] for v in row])
    image_minors = [image[i][a] * image[j][b] - image[i][b] * image[j][a]
                    for i, j in combinations(range(5), 2)
                    for a, b in combinations(range(5), 2)]
    index = reduce(gcd, (abs(v) for v in image_minors))
    assert index == s_gcd * agcd // (gval * prod(divisors))
    assert index == epsilon * delta * rho
    assert 3 * bproduct ** 3 % index == 0

    def circuit(support):
        vector = [0] * 5
        for k, i in enumerate(support):
            triple = tuple(j for j in support if j != i)
            vector[i] = (-1) ** k * aminors[triple] // L
        return vector
    first, second = circuit((0, 1, 2, 3)), circuit((0, 1, 2, 4))
    ell = next([a + t * b for a, b in zip(first, second)]
               for t in range(1, 5) if all(a + t * b for a, b in zip(first, second)))
    target = [index * e for e in ell]
    row_pair = next((i, j) for i, j in combinations(range(5), 2)
                    if any(image[i][a] * image[j][b] != image[i][b] * image[j][a]
                           for a, b in combinations(range(5), 2)))
    solution = solve_two_rows([image[i] for i in row_pair], [target[i] for i in row_pair])
    assert [sum(a * b for a, b in zip(row, solution)) for row in image] == target
    q = [sum(kernel[j][i] * solution[j] for j in range(5)) for i in range(6)]
    assert sum(a * b for a, b in zip(evaluations, q)) == 0
    assert all(sum(a * b for a, b in zip(jets[i], q)) == divisors[i] * target[i]
               for i in range(5))
    assert all(target)
    return epsilon


def main():
    basis, product_basis, jet_coefficients = quadratic_product_data()
    epsilons = [check_source(basis, product_basis, jet_coefficients, exponents)
                for exponents in ((0, 0, 0, 0, 0), (1, 1, 1, 1, 1), (0, 1, 2, 1, 3))]
    print("PASS: 20 first-jet sections, rank16; integral quadratic-product basis and exact index3.")
    print("PASS: all 200 complementary-triangle identities as complete chart polynomials.")
    print("PASS: 600 literal complementary-triangle first-jet minors and three exact lift indices.")
    print("PASS: delta divides product(b)^2; rho divides product(b); all forced core factors cancel.")
    print("PASS: three explicit common-Q lifts of full-support cofactor stresses.")
    print("No coefficient-height bound on those lifts or uniform arc bound is asserted.")


if __name__ == "__main__":
    main()
