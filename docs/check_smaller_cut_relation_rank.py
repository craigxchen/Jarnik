"""Exact algebra for smaller_cut_relation_rank.md (standard library only)."""

from itertools import combinations
from math import comb


def add(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, 0) + coefficient
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def scale(poly, coefficient):
    return {monomial: coefficient * value for monomial, value in poly.items() if coefficient * value}


def mul(left, right):
    out = {}
    for a, ca in left.items():
        for b, cb in right.items():
            monomial = tuple(x + y for x, y in zip(a, b))
            out[monomial] = out.get(monomial, 0) + ca * cb
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def variable(n, j):
    return {tuple(int(k == j) for k in range(n)): 1}


def difference(n, i, j):
    return add(variable(n, i), scale(variable(n, j), -1))


def derivative(poly, j):
    out = {}
    for monomial, coefficient in poly.items():
        if monomial[j]:
            reduced = list(monomial)
            reduced[j] -= 1
            out[tuple(reduced)] = coefficient * monomial[j]
    return out


def coefficient(poly, n, i, j):
    return poly.get(tuple(int(k == i) + int(k == j) for k in range(n)), 0)


def check_quadratic_bases():
    for n in range(4, 13):
        pairs = [pair for pair in combinations(range(n - 1), 2) if pair != (0, 1)]
        assert len(pairs) == n * (n - 3) // 2
        reference = mul(difference(n, 0, n - 1), difference(n, 1, n - 1))
        for i, j in pairs:
            poly = add(mul(difference(n, i, n - 1), difference(n, j, n - 1)),
                       scale(reference, -1))
            assert all(sum(monomial) == 2 and max(monomial) <= 1 for monomial in poly)
            total_derivative = {}
            for k in range(n):
                total_derivative = add(total_derivative, derivative(poly, k))
            assert not total_derivative
            assert [coefficient(poly, n, a, b) for a, b in pairs] == [
                int((a, b) == (i, j)) for a, b in pairs
            ]
            # Each basis element is a sum of at most two four-distinct-row matchings.
            if i == 0:
                expansion = mul(difference(n, 0, n - 1), difference(n, j, 1))
            elif i == 1:
                expansion = mul(difference(n, 1, n - 1), difference(n, j, 0))
            else:
                expansion = add(mul(difference(n, i, 0), difference(n, j, n - 1)),
                                mul(difference(n, 0, n - 1), difference(n, j, 1)))
            assert poly == expansion
        # The distinguished element q_13 is one nonzero matching product.
        distinguished = add(mul(difference(n, 0, n - 1), difference(n, 2, n - 1)),
                            scale(reference, -1))
        assert distinguished == mul(difference(n, 0, n - 1), difference(n, 2, 1))


def bracket(n, i, j, inside):
    if i in inside and j in inside:
        return {}
    if i in inside:
        return {(0,) * n: 1}
    if j in inside:
        return {(0,) * n: -1}
    return difference(n, i, j)


def restrict_graph(edges, inside, n=10):
    poly = {(0,) * n: 1}
    for i, j in edges:
        poly = mul(poly, bracket(n, i, j, inside))
    return poly


def check_common_factor():
    triangles = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]
    q_edges = triangles + [(6, 7)] * 2 + [(8, 9)] * 2
    r_edges = triangles + [(6, 8)] * 2 + [(7, 9)] * 2
    nonzero = 0
    for cut in combinations(range(10), 4):
        first = restrict_graph(q_edges, set(cut))
        second = restrict_graph(r_edges, set(cut))
        monomials = sorted(set(first) | set(second))
        for a, b in combinations(monomials, 2):
            assert first.get(a, 0) * second.get(b, 0) == first.get(b, 0) * second.get(a, 0)
        nonzero += bool(first or second)
    assert nonzero > 0
    # The residual factors are independent: the first vanishes at z7=z8,
    # whereas the second is nonzero if z7=z8=0,z9=1,z10=2.
    residual_points = {6: 0, 7: 0, 8: 1, 9: 2}
    q_value = (residual_points[6] - residual_points[7])**2 * (residual_points[8] - residual_points[9])**2
    r_value = (residual_points[6] - residual_points[8])**2 * (residual_points[7] - residual_points[9])**2
    assert q_value == 0 and r_value != 0
    print("All 210 first-smaller cuts of the independent ten-row pair have restriction rank at most one.")


def invariant_dimension(m):
    coefficients = [1]
    for _ in range(m):
        updated = [0] * (len(coefficients) + 2)
        for j, value in enumerate(coefficients):
            for k in range(3):
                updated[j + k] += value
        coefficients = updated
    return coefficients[m] - coefficients[m - 1]


def check_threshold_table():
    print("m, dim V, A_m, local dimension t, t*A_m/(dim V-t):")
    for m in (18, 24, 30, 32, 36):
        q = m // 2
        r = invariant_dimension(m)
        height = m * 2**(m - 2) - q * comb(m, q)
        t = (q + 1) * (q - 2) // 2
        print(m, r, height, t, f"{t * height / (r - t):.6f}")
    m = 30
    q = m // 2
    r = invariant_dimension(m)
    height = m * 2**(m - 2) - q * comb(m, q)
    t = (q + 1) * (q - 2) // 2
    balanced = comb(m, q) // 2
    cuts = comb(m, q - 1)
    assert t * height < 2 * (r - t)
    assert 2 * (balanced + cuts * t) < height * t


def main():
    check_quadratic_bases()
    print("Integral quadratic bases, coefficient projections, and four-cycle decompositions pass for n=4,...,12.")
    check_common_factor()
    check_threshold_table()
    print("No numerical zero relation or actual full-profile lattice construction is asserted.")


if __name__ == "__main__":
    main()
