"""A fixed-Q clean-prime content obstruction, with correction growth.

Only one core grows. This is not a full-fair or endpoint family.
"""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, prod

from check_five_row_first_jet_lifting import chart
from check_five_row_gradient_cofactor_stress import (
    bracket, conj, gaussian_gcd, mul, norm, power,
)
from check_five_row_gradient_discriminant_core_count import (
    evaluate, independent_graphs, quadratic_coefficients,
)
from check_five_row_petersen_circulation import EDGES, PAIRS, coefficient_table
from check_smaller_cut_relation_rank import add, derivative, scale


COEFFICIENTS = (35, 50, -4, -4, -4, -4)
PRIME = 29


def trim(poly):
    while poly and not poly[-1]:
        poly.pop()
    return poly


def remainder(first, second):
    first = trim(list(map(Fraction, first)))
    second = trim(list(map(Fraction, second)))
    while first and len(first) >= len(second):
        coefficient = first[-1] / second[-1]
        offset = len(first) - len(second)
        for j, c in enumerate(second):
            first[j + offset] -= coefficient * c
        trim(first)
    return first


def polynomial_gcd(first, second):
    first = trim(list(map(Fraction, first)))
    second = trim(list(map(Fraction, second)))
    while second:
        first, second = second, remainder(first, second)
    return first


def polynomial_gcd_mod(first, second, prime):
    first = trim([c % prime for c in first])
    second = trim([c % prime for c in second])
    while second:
        current = first[:]
        while current and len(current) >= len(second):
            coefficient = current[-1] * pow(second[-1], -1, prime) % prime
            offset = len(current) - len(second)
            for j, c in enumerate(second):
                current[j + offset] = (current[j + offset] - coefficient * c) % prime
            trim(current)
        first, second = second, current
    return first


def divide_gaussian(first, second):
    numerator = mul(first, conj(second))
    denominator = norm(second)
    assert all(c % denominator == 0 for c in numerator)
    return tuple(c // denominator for c in numerator)


def valuation(integer, prime):
    integer = abs(integer)
    assert integer
    count = 0
    while integer % prime == 0:
        count += 1
        integer //= prime
    return count


def fixed_section():
    basis = independent_graphs(5, 2)
    q = {}
    for coefficient, (_, poly) in zip(COEFFICIENTS, basis):
        q = add(q, scale(poly, coefficient))
    assert sum(map(abs, q.values())) == 1884
    chart_polynomial = chart(q)
    assert chart_polynomial == {
        (1, 2): 19, (2, 1): -27, (2, 0): 81, (0, 2): 15,
        (0, 1): 35, (1, 0): -39, (1, 1): -84,
    }
    A, B, C = [trim([chart_polynomial.get((a_degree, b_degree), 0)
                     for b_degree in range(3)]) for a_degree in (2, 1, 0)]
    assert (A, B, C) == ([81, -27], [-39, -84, 19], [0, 35, 15])
    convolution = lambda p, q: [sum(p[i] * q[k - i]
                                    for i in range(len(p)) if 0 <= k - i < len(q))
                               for k in range(len(p) + len(q) - 1)]
    BB, AC = convolution(B, B), convolution(A, C)
    discriminant = [BB[k] - 4 * (AC[k] if k < len(AC) else 0) for k in range(5)]
    assert discriminant == [1521, -4788, 4494, -1572, 361]
    discriminant_derivative = [i * discriminant[i] for i in range(1, 5)]
    assert len(polynomial_gcd(discriminant, discriminant_derivative)) == 1
    assert discriminant[-1] % PRIME != 0
    assert len(polynomial_gcd_mod(discriminant, discriminant_derivative, PRIME)) == 1
    assert len(polynomial_gcd(polynomial_gcd(A, B), C)) == 1
    table = coefficient_table([edges for edges, _ in basis])
    flows = []
    for i, j in EDGES:
        first, second = PAIRS[i], PAIRS[j]
        remaining = next(k for k in range(5) if k not in first + second)
        flows.append(sum(COEFFICIENTS[k] * table[first][k][remaining]
                         for k in range(6)))
    assert flows == [8, -4, -4, -35, 39, -4, 27, -81, 54, -19, -15, 34, 8, -50, -42]
    assert all(flows)
    assert all(value % PRIME for value in flows)
    rows = ((0, 1), (1, 0), (1, 1), (1, 2), (1, 3))
    values = sum(rows, ())
    assert evaluate(q, values) == 0
    multipliers = []
    for i, (x, y) in enumerate(rows):
        gx = evaluate(derivative(q, 2 * i), values)
        gy = evaluate(derivative(q, 2 * i + 1), values)
        lam = gy // x if x else -gx // y
        assert (gx, gy) == (-lam * y, lam * x)
        multipliers.append(lam)
    assert multipliers == [-222, 34, 9, -120, 77]
    assert all(value % PRIME for value in multipliers)
    assert gcd(120, PRIME) == 1
    return q


def prime_power_fixture(q, depth):
    H = power((5, 2), depth)
    n = norm(H)
    assert n == PRIME ** depth
    residue = (-H[0] * pow(H[1], -1, n)) % n
    choices = []
    for k in range(2 * PRIME):
        r = residue + k * n
        if r % 2:
            continue
        rows = [(-r - 2 * n * t, 1) for t in range(4)]
        corrections = [divide_gaussian(p, H) for p in rows]
        if all(norm(c) % PRIME for c in corrections):
            choices.append((r, rows, corrections))
    assert choices
    r, inside_rows, corrections = choices[0]
    assert 0 <= r < 58 * n
    rows = [(-1, 0)] + inside_rows
    assert all(norm(gaussian_gcd(p, conj(p))) == 1 for p in rows)
    values = sum(rows, ())
    assert evaluate(q, values) == 0
    for i, (x, y) in enumerate(rows):
        gx = evaluate(derivative(q, 2 * i), values)
        gy = evaluate(derivative(q, 2 * i + 1), values)
        assert gx * x + gy * y == 0 and (gx or gy)
    b_factors = []
    residues = []
    for i, j in combinations(range(5), 2):
        determinant = bracket(rows[i], rows[j])
        assert determinant
        core_norm = n if i and j else 1
        actual_gcd_norm = norm(gaussian_gcd(rows[i], rows[j]))
        assert actual_gcd_norm == core_norm
        assert determinant % actual_gcd_norm == 0
        residues.append(determinant // actual_gcd_norm)
        b_factors.append(abs(determinant) // core_norm)
    assert max(map(abs, residues)) <= 6
    assert prod(b_factors) == 768
    assert all(b % PRIME for b in b_factors)

    cpoly, bpoly, apoly = quadratic_coefficients(q, 3)
    others = values[:6] + values[8:]
    a, b, c = [evaluate(poly, others) for poly in (apoly, bpoly, cpoly)]
    assert (a, b, c) == (0, 480 * n * n, 480 * n * n * (r + 4 * n))
    content = reduce(gcd, (abs(a), abs(b), abs(c)))
    assert content == 480 * n * n
    F_i = n
    K_i_norm = norm(corrections[2])
    assert K_i_norm % PRIME != 0
    excess = content * K_i_norm // F_i
    assert excess == 480 * n * K_i_norm
    assert valuation(content, PRIME) == 2 * depth
    assert valuation(excess, PRIME) == depth
    # Uniform archimedean bounds make the growing correction cost explicit.
    assert 36 * n <= max(map(norm, corrections)) <= 4097 * n
    assert 16 * n <= K_i_norm <= 3845 * n
    assert 7680 * n * n <= excess <= 1845600 * n * n
    assert gcd(120, PRIME) == 1  # normalized coefficient content at b=3


if __name__ == "__main__":
    q = fixed_section()
    for depth in range(1, 7):
        prime_power_fixture(q, depth)
    print("PASS: fixed CQ=1884 section; discriminant derived from chart and squarefree over Q and mod29.")
    print("PASS: all15 boundary-node values, all5 normalized old gradients, and gcd120 are29-units.")
    print("PASS: six clean29-power frames, five unramified old rows, all actual residues bounded6.")
    print("PASS: exact g=480n^2, F=n, v29(E)=m, and explicit growing correction bounds.")
    print("This is not a full-fair family or a counterexample to a bound charging exp(B sigma).")
