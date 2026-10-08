#!/usr/bin/env python3
"""Exact algebra/interval checks for the positive seven-radical defect."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod

from check_gale_square_quartic_characterization import primitive_pluckers


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    while len(out) > 1 and not out[-1]:
        out.pop()
    return out


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return add(out, [])


def derivative(a):
    return [i * a[i] for i in range(1, len(a))] or [F(0)]


def check_sixth_derivative(form):
    a, b, c = map(F, form)
    q = [a, 2 * b, c]
    current = [F(1)]
    for k in range(6):
        current = add(multiply(q, derivative(current)),
                      [v * (F(5, 2) - k) for v in multiply(derivative(q), current)])
    assert current == [225 * (a * c - b * b) ** 3]


def rational_gcd(values):
    denominator = lcm(*(value.denominator for value in values))
    return F(reduce(gcd, (int(value * denominator) for value in values)), denominator)


def sqrt_interval(integer, scale=10**70):
    lower = isqrt(integer * scale * scale)
    if lower * lower == integer * scale * scale:
        return F(lower, scale), F(lower, scale)
    return F(lower, scale), F(lower + 1, scale)


def linear_interval(coefficients, radicands):
    lower = upper = F(0)
    for coefficient, radicand in zip(coefficients, radicands):
        lo, hi = sqrt_interval(radicand)
        lower += coefficient * (lo if coefficient >= 0 else hi)
        upper += coefficient * (hi if coefficient >= 0 else lo)
    return lower, upper


def bezout(values):
    def extended(a, b):
        if not b:
            return a, 1, 0
        g, x, y = extended(b, a % b)
        return g, y, x - (a // b) * y
    g = 0
    coefficients = []
    for value in values:
        g, x, y = extended(g, abs(value))
        coefficients = [x * coefficient for coefficient in coefficients]
        coefficients.append(y * (1 if value >= 0 else -1))
    assert sum(a*b for a, b in zip(coefficients, values)) == g
    return g, coefficients


def check_saturated_basis_scaling(pluckers, cs, star_content):
    def bracket(i, j):
        return 0 if i == j else (pluckers[i, j] if i < j else -pluckers[j, i])
    first_content = reduce(gcd, (abs(bracket(0, j)) for j in range(1, 7)))
    ys = [bracket(0, j) // first_content for j in range(7)]
    g, coefficients = bezout(ys)
    assert g == 1
    xs = [sum(coefficients[k] * bracket(j, k) for k in range(7)) for j in range(7)]
    shift = next(k for k in range(8) if all(x + k*y for x, y in zip(xs, ys)))
    basis = [(x + shift*y, y) for x, y in zip(xs, ys)]
    assert all(sum(row[k] for row in basis) == 0 for k in (0, 1))
    for i, j in combinations(range(7), 2):
        assert basis[i][0]*basis[j][1] - basis[i][1]*basis[j][0] == bracket(i, j)
    # Primitive minors certify that these integral columns are saturated.
    assert reduce(gcd, map(abs, pluckers.values())) == 1
    matrix = [[F(x*x), F(2*x*y), F(y*y), F(cs[i])]
              for i, (x, y) in enumerate(basis[:3])]
    for column in range(3):
        pivot = next(i for i in range(column, 3) if matrix[i][column])
        matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        scale = matrix[column][column]
        matrix[column] = [x / scale for x in matrix[column]]
        for i in range(3):
            if i != column:
                scale = matrix[i][column]
                matrix[i] = [x - scale*y for x, y in zip(matrix[i], matrix[column])]
    a, b, c = [matrix[i][3] for i in range(3)]
    assert a > 0 and a*c-b*b > 0
    assert all(a*x*x + 2*b*x*y + c*y*y == value for (x, y), value in zip(basis, cs))
    contents = [gcd(abs(x), abs(y)) for x, y in basis]
    signs = [1 if x > 0 else -1 for x, _ in basis]
    marked = [(sign*x//ell, sign*y//ell) for (x, y), sign, ell in zip(basis, signs, contents)]
    values = [a*x*x + 2*b*x*y + c*y*y for x, y in marked]
    determinants = [prod(v[0]*w[1]-v[1]*w[0] for j, w in enumerate(marked) if i != j)
                    for i, v in enumerate(marked)]
    row_product = prod(contents)
    kappa = F(-prod(signs) * star_content**2, row_product)
    assert all(kappa * q*q / d == sign*ell
               for q, d, sign, ell in zip(values, determinants, signs, contents))
    tau = rational_gcd([q**5 / d**2 for q, d in zip(values, determinants)])
    assert tau == F(row_product**2, star_content**4) == 1/kappa**2
    assert all((1 if d > 0 else -1) == -prod(signs)*sign for d, sign in zip(determinants, signs))
    return basis, (a, b, c), contents, signs


def check_configuration(form, vectors, all_squares=False):
    a, b, c = map(F, form)
    vectors = [tuple(-z for z in v) if v[0] < 0 else v for v in vectors]
    assert all(x > 0 for x, _ in vectors)
    det = lambda v, w: v[0] * w[1] - v[1] * w[0]
    values = [a*x*x + 2*b*x*y + c*y*y for x, y in vectors]
    stars = [prod(det(v, w) for j, w in enumerate(vectors) if i != j)
             for i, v in enumerate(vectors)]
    rho_squared = [q**5 / d**2 for q, d in zip(values, stars)]
    tau = rational_gcd(rho_squared)
    cs = [int(value / tau) for value in rho_squared]
    assert all(F(c) == value / tau for c, value in zip(cs, rho_squared))
    assert reduce(gcd, cs) == 1
    signs = [1 if d > 0 else -1 for d in stars]
    gale = [[q**2 / d * z for z in v] for q, d, v in zip(values, stars, vectors)]
    assert all(sum(row[k] for row in gale) == 0 for k in (0, 1))
    _, pluckers = primitive_pluckers(gale)
    def bracket(i, j):
        return pluckers[i, j] if i < j else -pluckers[j, i]
    roots = [isqrt(-prod(bracket(i, j) for j in range(7) if i != j)) for i in range(7)]
    content = reduce(gcd, roots)
    assert cs == [root // content for root in roots]
    check_saturated_basis_scaling(pluckers, cs, content)
    parameters = [F(y, x) for x, y in vectors]
    left, right = min(parameters), max(parameters)
    evaluate = lambda t: a + 2*b*t + c*t*t
    candidates = [evaluate(left), evaluate(right)]
    if left <= -b/c <= right:
        candidates.append(evaluate(-b/c))
    q_min, q_max = min(candidates), max(candidates)
    x_product = prod(x for x, _ in vectors)
    squared_constant = F(25, 256) * (a*c-b*b)**6 / (x_product**2 * tau)
    lo, hi = linear_interval(signs, cs)
    assert lo > 0
    assert lo * lo >= squared_constant / q_max**7
    assert hi * hi <= squared_constant / q_min**7
    if all_squares:
        assert all(isqrt(c)**2 == c for c in cs)
        assert lo == hi and lo.denominator == 1 and lo >= 1
    return len(str(max(cs)))


def check_fixed_radical_pigeonhole():
    radicands = [1, 5, 13, 17, 65, 85, 221]
    scale = 10**70
    approximations = [isqrt(b * scale * scale) for b in radicands[1:]]
    height = 4
    values = []
    for coefficients in product(range(height + 1), repeat=6):
        value = sum(a*b for a, b in zip(coefficients, approximations))
        values.append((value % scale, coefficients))
    values.sort()
    index = min(range(len(values) - 1), key=lambda i: values[i + 1][0] - values[i][0])
    lower, upper = values[index], values[index + 1]
    coefficients = [b - a for a, b in zip(lower[1], upper[1])]
    raw = sum(a*b for a, b in zip(coefficients, approximations))
    constant = -(raw // scale)
    coefficients = [constant] + coefficients
    lo, hi = linear_interval(coefficients, radicands)
    assert 0 < lo <= hi < F(1, len(values) - 1)
    assert all(abs(a) <= height for a in coefficients[1:])
    print("PASS: certified fixed-radical small form", coefficients, "; upper bound 1/15624.")


if __name__ == "__main__":
    forms = [(1, 0, 1), (1, 2, 5), (2, 1, 5), (F(2, 3), F(1, 3), F(5, 3))]
    node_sets = [[(1, i) for i in range(-3, 4)], [(2, i) for i in range(-3, 4)]]
    for form in forms:
        check_sixth_derivative(form)
        for vectors in node_sets:
            check_configuration(form, vectors)
    square_vectors = [(u*u - 1, 2*u) for u in (0, 2, 3, 4, 5, 6, 7)]
    check_configuration((1, 0, 1), square_vectors, all_squares=True)
    print("PASS: exact sixth derivative; nine intrinsic-normalization and positive-bound configurations.")
    print("PASS: saturated integer basis reconstructions, orientations, zero sums, and exact K/tau scaling.")
    print("PASS: integer-radical case retains the essential rational gcd normalization.")
    check_fixed_radical_pigeonhole()
