"""Exact Vieta-swap arithmetic on literal five-row Gaussian sources.

The full-cut fixtures have unrestricted invariant height and are not
endpoint arcs.  They audit content, primitive roots, pair residues,
core transfer, and Gaussian anchor divisibility only.
"""

from math import gcd, prod

from check_five_row_gradient_cofactor_stress import (
    CUTS, bracket, conj, gaussian_gcd, mul, norm, popcount, power,
    split_primes,
)
from check_five_row_gradient_discriminant_core_count import independent_graphs
from check_smaller_cut_relation_rank import add as poly_add, mul as poly_mul, scale as poly_scale


def graph_value(edges, rows):
    return prod(bracket(rows[i], rows[j]) for i, j in edges)


def extended_gcd(a, b):
    if b == 0:
        return abs(a), (1 if a > 0 else -1), 0
    g, x, y = extended_gcd(b, a % b)
    return g, y, x - (a // b) * y


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(a, n):
    return n * a[0], n * a[1]


def quotient(a, b):
    numerator = mul(a, conj(b))
    denominator = norm(b)
    assert denominator and all(n % denominator == 0 for n in numerator)
    return numerator[0] // denominator, numerator[1] // denominator


def divides(a, b):
    try:
        quotient(b, a)
        return True
    except AssertionError:
        return False


def orders(mask, i):
    a = popcount(mask & ~(1 << i))
    inside = bool(mask & (1 << i))
    h = max(0, 2 * a - (5 if inside else 4))
    selected = (not inside and a >= 2) or (inside and a >= 3)
    return h, selected


def check_formal_counts():
    for i in range(5):
        assert sum(orders(mask, i)[0] for mask in CUTS) == 19
        selected = [mask for mask in CUTS if orders(mask, i)[1]]
        assert len(selected) == 16
        for j in range(5):
            if j != i:
                assert sum(bool(mask & (1 << j)) for mask in selected) == 11


def check_collision_restrictions_and_matching_counts():
    one, zero, t = {(0,): 1}, {}, {(1,): 1}
    # Collision 0=1, with retained directions (infinity,0,1,t).
    representatives = ((zero, one), (zero, one), (one, zero),
                       (one, one), (one, t))

    def boundary_bracket(i, j):
        xi, yi = representatives[i]
        xj, yj = representatives[j]
        return poly_add(poly_mul(xi, yj), poly_scale(poly_mul(yi, xj), -1))

    restricted = []
    for edges, _ in independent_graphs(5, 2):
        poly = one
        for i, j in edges:
            poly = poly_mul(poly, boundary_bracket(i, j))
        restricted.append(poly)
    assert restricted == [zero, zero, zero, zero,
                          {(0,): -1, (1,): 1}, {(1,): 1}]

    retained = (1, 2, 3, 4)
    matchings = (((1, 2), (3, 4)), ((1, 3), (2, 4)),
                 ((1, 4), (2, 3)))
    common_weight = 0
    for mask in CUTS:
        a = sum(bool(mask & (1 << j)) for j in retained)
        minimum = min(sum(bool(mask & (1 << i)) and bool(mask & (1 << j))
                          for i, j in matching)
                      for matching in matchings)
        assert minimum == max(0, a - 2)
        common_weight += minimum
    assert common_weight == 12
    for matching in matchings:
        assert sum(sum(bool(mask & (1 << i)) and bool(mask & (1 << j))
                       for i, j in matching)
                   for mask in CUTS) == 16


def literal_source(correction_exponents):
    primes = split_primes(32)
    blocks = dict(zip(CUTS, primes[:31]))
    weights = {mask: norm(blocks[mask]) for mask in CUTS}
    corrections = [power(primes[31], exponent)
                   for exponent in correction_exponents]
    rows = []
    for i in range(5):
        p = corrections[i]
        for mask in CUTS:
            if mask & (1 << i):
                p = mul(p, blocks[mask])
        assert gcd(*p) == 1
        assert norm(gaussian_gcd(p, conj(p))) == 1
        rows.append(p)

    # A literal numerical zero.  The graph coefficients are large in
    # these fixtures, so this is not a small-invariant endpoint test.
    graph_basis = independent_graphs(5, 2)
    first, second = graph_basis[0][0], graph_basis[1][0]
    first_value = graph_value(first, rows)
    second_value = graph_value(second, rows)
    divisor = gcd(first_value, second_value)
    coefficients = second_value // divisor, -first_value // divisor
    assert coefficients != (0, 0)

    def q_with_row(i, replacement):
        changed = list(rows)
        changed[i] = replacement
        return (coefficients[0] * graph_value(first, changed)
                + coefficients[1] * graph_value(second, changed))

    assert all(q_with_row(i, rows[i]) == 0 for i in range(5))
    source_anchor = (1, 0)
    for p in rows:
        source_anchor = mul(source_anchor, conj(p))

    for i, p in enumerate(rows):
        a = q_with_row(i, (1, 0))
        c = q_with_row(i, (0, 1))
        b = q_with_row(i, (1, 1)) - a - c
        content = gcd(a, gcd(b, c))
        assert content > 0
        x, y = p
        g, bezout_x, bezout_y = extended_gcd(x, y)
        assert g == 1
        frame = (-bezout_y, bezout_x)
        assert bracket(p, frame) == 1
        lam = (2 * a * x + b * y) * frame[0] + (b * x + 2 * c * y) * frame[1]
        assert lam != 0
        C = q_with_row(i, frame)
        assert content == gcd(lam, C)
        assert b * b - 4 * a * c == lam * lam
        # Scale by 1/content using exact Cartesian divisibility.
        numerator = add(scale(p, C), scale(frame, -lam))
        assert numerator[0] % content == numerator[1] % content == 0
        swapped = (numerator[0] // content, numerator[1] // content)
        assert gcd(*swapped) == 1 and bracket(p, swapped) == -lam // content
        assert q_with_row(i, swapped) == 0
        assert (c - a, -b) == scale(mul(p, swapped), content)
        for shift in (-3, 4):
            shifted_frame = add(frame, scale(p, shift))
            shifted_C = q_with_row(i, shifted_frame)
            assert shifted_C == C + shift * lam
            assert add(scale(p, shifted_C), scale(shifted_frame, -lam)) == numerator

        F = prod(weights[mask] ** orders(mask, i)[0] for mask in CUTS)
        assert content * norm(corrections[i]) % F == 0
        E = content * norm(corrections[i]) // F
        assert E > 0
        selected = [mask for mask in CUTS if orders(mask, i)[1]]
        transfer = (1, 0)
        for mask in selected:
            transfer = mul(transfer, blocks[mask])
        assert divides(transfer, scale(swapped, E))
        lost = quotient(transfer, gaussian_gcd(transfer, swapped))
        assert E % norm(lost) == 0
        assert divides(conj(transfer), source_anchor)

        for j, other in enumerate(rows):
            if j == i:
                continue
            old_delta = bracket(p, other)
            new_delta = bracket(swapped, other)
            assert old_delta != 0
            assert q_with_row(i, other) == content * old_delta * new_delta
            common = norm(gaussian_gcd(swapped, other))
            assert new_delta % common == 0
            shared_core = prod(weights[mask] for mask in selected
                               if mask & (1 << j))
            assert common * E >= shared_core

        reduced_gcd = gaussian_gcd(swapped, conj(swapped))
        assert norm(reduced_gcd) in (1, 2)
        numerator_phase = quotient(swapped, reduced_gcd)
        denominator_phase = quotient(conj(swapped), reduced_gcd)
        assert norm(gaussian_gcd(numerator_phase, denominator_phase)) == 1
        anchor_multiplier = quotient(
            denominator_phase,
            gaussian_gcd(denominator_phase, source_anchor),
        )
        assert (norm(gaussian_gcd(denominator_phase, source_anchor)) * E
                >= norm(transfer))
        assert norm(anchor_multiplier) * norm(transfer) <= norm(swapped) * E
        assert divides(denominator_phase,
                       mul(anchor_multiplier, source_anchor))
        new_anchor = mul(anchor_multiplier, source_anchor)
        assert divides(denominator_phase, new_anchor)
        new_point = mul(quotient(new_anchor, denominator_phase), numerator_phase)
        assert mul(new_point, denominator_phase) == mul(new_anchor, numerator_phase)


def small_algebraic_fixture():
    rows = [(1, 1), (2, 1), (1, 2), (3, 2), (2, 3)]
    basis = independent_graphs(5, 2)
    q_coefficients = (-15, -4)
    polynomial = poly_add(poly_scale(basis[0][1], q_coefficients[0]),
                          poly_scale(basis[1][1], q_coefficients[1]))
    assert sum(abs(coefficient) for coefficient in polynomial.values()) == 440

    def q_with_row(i, replacement):
        changed = list(rows)
        changed[i] = replacement
        return sum(q_coefficients[k] * graph_value(basis[k][0], changed)
                   for k in (0, 1))

    assert q_with_row(0, rows[0]) == 0
    assert q_with_row(2, rows[2]) == 0
    # The canonical swap at row two changes squared norm 5 to 149;
    # the canonical swap at row zero collides projectively with row one.
    for i, expected_norm in ((0, 5), (2, 149)):
        p = rows[i]
        a = q_with_row(i, (1, 0))
        c = q_with_row(i, (0, 1))
        b = q_with_row(i, (1, 1)) - a - c
        g, sx, sy = extended_gcd(*p)
        assert g == 1
        frame = (-sy, sx)
        lam = (2 * a * p[0] + b * p[1]) * frame[0] + (b * p[0] + 2 * c * p[1]) * frame[1]
        C = q_with_row(i, frame)
        content = gcd(lam, C)
        root = add(scale(p, C // content), scale(frame, -lam // content))
        assert norm(root) == expected_norm
        if i == 0:
            assert bracket(root, rows[1]) == 0


if __name__ == "__main__":
    check_formal_counts()
    check_collision_restrictions_and_matching_counts()
    literal_source((0, 0, 0, 0, 0))
    literal_source((0, 1, 2, 1, 3))
    small_algebraic_fixture()
    print("PASS: 19-weight content and 16/11 cut transfer counts; ten literal full-cut Vieta swaps.")
    print("PASS: primitive second roots, exact pair residues, correction E cancellation, and Gaussian anchor divisibility.")
    print("PASS: nonzero collision boundary restriction basis and four-row matching core weights 16/12.")
    print("PASS: low-height algebraic norm-change and collision fixtures; no endpoint or descent claim.")
