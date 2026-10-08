"""Exact cut and rational-fixture checks for the five-swap cross-ratio route.

The universal sum-six identity is proved by the companion symbolic
certificate.  These fixtures exercise it on literal integer invariant
sections and audit the formal core valuations; they are not endpoint
sources and impose no small residue budget.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, lcm

from check_five_row_gradient_discriminant_core_count import (
    evaluate, independent_graphs, inverse,
)
from check_five_row_simultaneous_vieta_identity import affine
from check_smaller_cut_relation_rank import derivative, restrict_graph


def popcount(mask):
    return bin(mask).count("1")


def cut_term_exponents():
    results = {}
    for mask in range(1, 32):
        s = popcount(mask)
        old = [int(bool(mask & (1 << i))) for i in range(5)]
        selected = [int((not old[i] and s >= 2) or (old[i] and s >= 4))
                    for i in range(5)]
        self_exponent = [int(old[i] and s >= 4) for i in range(5)]
        exponents = []
        for i in range(5):
            for j in range(i + 1, 5):
                exponent = (old[i] * old[j] + selected[i] * selected[j]
                            - self_exponent[i] - self_exponent[j])
                exponents.append(exponent)
        expected = {1: (10, 0), 2: (6, 4), 3: (6, 4),
                    4: (10, 0), 5: (10, 0)}[s]
        assert (exponents.count(0), exponents.count(1)) == expected
        assert min(exponents) == 0
        results[mask] = tuple(exponents)
    return results


def graph_value(edges, x):
    value = Fraction(1)
    for i, j in edges:
        value *= x[j] - x[i]
    return value


def integer_section_zero_at(x, edges):
    graph_values = [graph_value(graph, x) for graph in edges]
    trial = sum((j + 1) * graph_values[j] for j in range(1, 6))
    coefficients = [-trial] + [(j + 1) * graph_values[0]
                               for j in range(1, 6)]
    denominator = 1
    for value in coefficients:
        denominator = lcm(denominator, value.denominator)
    integers = [int(value * denominator) for value in coefficients]
    assert any(integers)
    assert sum(integers[j] * graph_values[j] for j in range(6)) == 0
    return integers


def canonical_binary_factor(a_X2, b_XY, c_Y2, old_y):
    """Check the signed Vieta factorization for P=(1,old_y)."""
    a_X2, b_XY, c_Y2, old_y = map(int, (a_X2, b_XY, c_Y2, old_y))
    assert a_X2 + b_XY * old_y + c_Y2 * old_y**2 == 0
    content = gcd(a_X2, b_XY, c_Y2)
    assert content > 0
    lam = b_XY + 2 * c_Y2 * old_y
    C = c_Y2
    assert content == gcd(lam, C)
    # u=(0,1) has determinant one with P=(1,old_y).
    new_X = C // content
    new_Y = (C * old_y - lam) // content
    assert gcd(new_X, new_Y) == 1
    assert a_X2 == content * old_y * new_Y
    assert b_XY == -content * (new_Y + old_y * new_X)
    assert c_Y2 == content * new_X
    assert new_Y - old_y * new_X == -lam // content
    if new_Y == 0:
        assert abs(new_X) == 1  # Phase P'/bar(P')=1: anchor collision.
    return new_X, new_Y


def swaps_and_crossratio(x, edges):
    coefficients = integer_section_zero_at(x, edges)

    def q_at_row(i, value):
        changed = list(x)
        changed[i] = value
        return sum(coefficients[j] * graph_value(edges[j], changed)
                   for j in range(6))

    swapped = []
    for i, old in enumerate(x):
        c = q_at_row(i, Fraction(0))
        a = (q_at_row(i, Fraction(2)) -
             2 * q_at_row(i, Fraction(1)) + c) / 2
        b = q_at_row(i, Fraction(1)) - a - c
        assert a and q_at_row(i, old) == 0
        new = -b / a - old
        assert new != old and q_at_row(i, new) == 0
        if all(value.denominator == 1 for value in x):
            new_X, new_Y = canonical_binary_factor(c, b, a, old)
            assert new == Fraction(new_Y, new_X)
        swapped.append(new)

    terms = [((x[i] - x[j]) * (swapped[i] - swapped[j]) /
              ((x[i] - swapped[i]) * (x[j] - swapped[j])))
             for i in range(5) for j in range(i + 1, 5)]
    assert sum(terms) == 6
    diameter = max(abs(x[i] - x[j])
                   for i in range(5) for j in range(i + 1, 5))
    assert min(abs(x[i] - swapped[i]) for i in range(5)) < 4 * diameter
    return tuple(swapped)


def one_near_only_same_section(T, basis):
    """Reconstruct a rational common Q for four escaping Vieta roots."""
    T = Fraction(T)
    x = tuple(map(Fraction, (0, 1, 2, 3, 4)))
    t = (-10 + 30 / T) / (6 - 10 / T - 20 / T**2)
    shifts = (t,) + (T,) * 4
    new = tuple(xi - ki for xi, ki in zip(x, shifts))
    terms = [((x[i] - x[j]) * (new[i] - new[j]) /
              (shifts[i] * shifts[j]))
             for i in range(5) for j in range(i + 1, 5)]
    assert sum(terms) == 6
    assert abs(t) < 16 and all(abs(k) > 16 for k in shifts[1:])

    u = tuple(1 / k for k in shifts)
    U = [sum(xi**r * ui for xi, ui in zip(x, u)) for r in range(3)]
    S = sum(x)
    assert (U[1] - 3) * (U[1] - 2) == U[0] * (U[2] - S)
    # A nonzero moment-kernel generator, and the required open condition.
    A, B = U[0], 3 - U[1]
    assert A or B
    assert all(A * xi + B for xi in x)
    d = [Fraction(1) for _ in x]
    for i in range(5):
        for j in range(5):
            if j != i:
                d[i] *= x[i] - x[j]
    lam = [(A * xi + B) / di for xi, di in zip(x, d)]
    a = [ui * li for ui, li in zip(u, lam)]
    matrix = [[evaluate(q, x) for q in basis]] + [
        [evaluate(derivative(derivative(q, i), i), x) / 2
         for q in basis] for i in range(5)]
    inv = inverse(matrix)
    targets = [Fraction(0)] + a
    coefficients = [sum(c * target for c, target in zip(row, targets))
                    for row in inv]
    denominator = 1
    for value in coefficients:
        denominator = lcm(denominator, value.denominator)
    integers = [int(value * denominator) for value in coefficients]
    content = gcd(*integers)
    assert content > 0
    integers = [value // content for value in integers]
    q = {}
    for coefficient, graph in zip(integers, basis):
        for monomial, value in graph.items():
            q[monomial] = q.get(monomial, 0) + coefficient * value
    q = {monomial: value for monomial, value in q.items() if value}
    assert q and evaluate(q, x) == 0
    for i in range(5):
        actual_lam = evaluate(derivative(q, i), x)
        actual_a = evaluate(derivative(derivative(q, i), i), x) / 2
        assert actual_lam and actual_a
        assert actual_lam / actual_a == shifts[i]
        changed = list(x)
        changed[i] = new[i]
        assert evaluate(q, tuple(changed)) == 0
    return t, q, integers


def mod11_squarefree(polynomial):
    def trim(values):
        values = [value % 11 for value in values]
        while values and values[-1] == 0:
            values.pop()
        return values

    def remainder(dividend, divisor):
        dividend, divisor = trim(dividend), trim(divisor)
        while len(dividend) >= len(divisor) and dividend:
            factor = dividend[-1] * pow(divisor[-1], -1, 11) % 11
            offset = len(dividend) - len(divisor)
            for j, value in enumerate(divisor):
                dividend[offset + j] = (dividend[offset + j] - factor * value) % 11
            dividend = trim(dividend)
        return dividend

    a = trim(polynomial)
    b = trim([j * polynomial[j] for j in range(1, len(polynomial))])
    assert len(a) == len(polynomial)  # Degree is preserved modulo 11.
    while b:
        a, b = b, remainder(a, b)
    return len(a) == 1


def irreducible_nonnodal_fixture(q, coefficients, edges):
    """Audit T=100 on the normalized M_0,5 affine chart."""
    assert coefficients == [-19948, 3804, -2018, -1843, 4038, -59]
    # Fix the first three projective rows at 0,1,2 and vary the last two.
    F = {}
    for monomial, coefficient in q.items():
        key = (monomial[3], monomial[4])
        F[key] = F.get(key, 0) + coefficient * 0**monomial[0] * 1**monomial[1] * 2**monomial[2]
    F = {key: value for key, value in F.items() if value}
    A = [F.get((j, 2), 0) for j in range(3)]
    B = [F.get((j, 1), 0) for j in range(3)]
    C = [F.get((j, 0), 0) for j in range(3)]
    assert A == [35388, -5871, -1961]
    assert B == [-71712, -18906, 15558]
    assert C == [0, 79320, -31816]

    def multiply(a, b):
        result = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            for j, bj in enumerate(b):
                result[i + j] += ai * bj
        return result

    b2, ac = multiply(B, B), multiply(A, C)
    disc = [b2[j] - 4 * ac[j] for j in range(5)]
    assert disc == [5142610944, -8516330496, 4492415556,
                    -713259960, -7513340]
    assert [value % 11 for value in disc] == [9, 8, 2, 5, 1]
    assert mod11_squarefree(disc)
    # The squarefree quartic discriminant is nonsquare in C(x), and
    # excludes a common nonconstant x-factor from A,B,C. Thus F is
    # irreducible over C[x,y] on the normalized interior chart.
    for pair in combinations(range(5), 2):
        restriction = {}
        for coefficient, graph in zip(coefficients, edges):
            for monomial, value in restrict_graph(graph, set(pair), 5).items():
                restriction[monomial] = (restriction.get(monomial, 0)
                                         + coefficient * value)
        outside = [restriction.get(tuple(int(i == j) for i in range(5)), 0)
                   for j in range(5) if j not in pair]
        assert sum(outside) == 0 and all(outside)
    # Each boundary restriction is a genuine linear three-term form;
    # no boundary component or boundary-node contact is present.


def main():
    assert len(cut_term_exponents()) == 31
    basis = [edges for edges, _ in independent_graphs(5, 2)]
    fixtures = (
        (0, 1, 2, 4, 7),
        (-3, -1, 0, 2, 5),
        (1, 2, 5, 9, 13),
        (Fraction(0), Fraction(1, 100), Fraction(2, 100),
         Fraction(4, 100), Fraction(7, 100)),
    )
    for fixture in fixtures:
        x = tuple(Fraction(value) for value in fixture)
        swaps_and_crossratio(x, basis)
    assert canonical_binary_factor(0, -3, 1, 3) == (1, 0)
    affine_basis = [affine(q) for _, q in independent_graphs(5, 2)]
    sharp = [one_near_only_same_section(T, affine_basis)
             for T in (100, 1000, 10000)]
    assert all(-2 < t < -1 for t, _, _ in sharp)
    irreducible_nonnodal_fixture(sharp[0][1], sharp[0][2], basis)
    print("PASS: thirty-one exact formal cut exponent vectors; six or ten terms have no forced core factor.")
    print("PASS: four rational old/new five-swap fixtures, sum-six identity and a near swap.")
    print("PASS: three rational common-Q reconstructions with exactly one near swap and four escaping swaps.")
    print("PASS: T=100 normalized discriminant squarefree modulo 11 and all ten boundary restrictions nonnodal.")
    print("PASS: signed X²-content factorization on integer sections, including a real-axis anchor collision.")
    print("No coefficient-height, anchor, or endpoint family claim is tested.")


if __name__ == "__main__":
    main()
