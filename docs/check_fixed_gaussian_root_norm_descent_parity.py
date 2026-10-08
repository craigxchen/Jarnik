"""Exact arithmetic checks for fixed-root norm descent and pair factoring.

The general field-degree parity obstruction is proved in the accompanying
note. No finite calculation here determines the coefficient field of the
certified four-row profile or proves an odd-degree exclusion.
"""

from fractions import Fraction
from itertools import combinations
from math import prod

from check_balanced_three_row_gaussian_polynomial_family import (
    FACTORS, INCIDENT, conjugate, mul, polynomial_product,
)


def trim(p):
    while len(p) > 1 and p[-1] == (0, 0):
        p.pop()
    return p


def add(p, q, scale=1):
    out = [(0, 0)] * max(len(p), len(q))
    for source, scalar in ((p, 1), (q, scale)):
        for j, (x, y) in enumerate(source):
            a, b = out[j]
            out[j] = a + scalar * x, b + scalar * y
    return trim(out)


def scale(p, k):
    return trim([(k * a, k * b) for a, b in p])


def conj(p):
    return [conjugate(z) for z in p]


def imag(p):
    return trim([(b, 0) for _, b in p])


def at_i(p):
    value = (0, 0)
    for a, b in reversed(p):
        x, y = mul(value, (0, 1))
        value = x + a, y + b
    return value


def derivative(p):
    return trim([(j * a, j * b) for j, (a, b) in enumerate(p)][1:])


def rational_remainder(p, q):
    assert all(b == 0 for _, b in p + q) and q != [(0, 0)]
    p = p[:]
    while p != [(0, 0)] and len(p) >= len(q):
        coefficient = Fraction(p[-1][0]) / q[-1][0]
        shifted = [(0, 0)] * (len(p) - len(q)) + scale(q, coefficient)
        p = add(p, shifted, -1)
    return p


def rational_gcd(p, q):
    while q != [(0, 0)]:
        p, q = q, rational_remainder(p, q)
    return scale(p, 1 / Fraction(p[-1][0]))


def product(polynomials):
    result = [(1, 0)]
    for p in polynomials:
        result = polynomial_product(result, p)
    return trim(result)


def real(coefficients):
    return [(x, 0) for x in coefficients]


def fixed_root_real_part(z, y):
    # X=(t^2+1)Z-Yt satisfies X(i)=-iY exactly.
    return add(polynomial_product(real([1, 0, 1]), z), real([0, -y]))


def check_parity():
    power = (1, 0)
    for d in range(1, 21):
        power = mul(power, (0, -2))
        # -(-2i)^d/(2i)=i(-2i)^d/2; integral for d>=1.
        forced = -power[1] // 2, power[0] // 2
        if d % 2:
            assert forced == ((-4) ** ((d - 1) // 2), 0)
        else:
            assert forced[0] == 0 and forced[1] != 0

    # All conjugate factors share i. These arbitrary rational factors
    # exercise the evaluation identity, without pretending they form a
    # field-conjugacy orbit or have a constant-imaginary product.
    for d in range(1, 12):
        ys = [j + 1 for j in range(d)]
        factors = []
        for j, y in enumerate(ys):
            x = fixed_root_real_part(real([j - 2, 1, 1]), y)
            factors.append(add(x, [(0, y)]))
        q = product(factors)
        assert at_i(q) == (0, 0)
        expected = (prod(ys), 0)
        for _ in range(d):
            expected = mul(expected, (0, -2))
        assert at_i(conj(q)) == expected
        if d % 2:
            residue = real([(-4) ** ((d - 1) // 2) * prod(ys)])
        else:
            residue = real([0, -2 * (-4) ** ((d - 2) // 2) * prod(ys)])
        assert rational_remainder(imag(q), real([1, 0, 1])) == residue


def check_quadratic():
    count = 0
    for c in range(-2, 3):
        for e in range(-2, 3):
            if c == e == 0:
                continue
            for k in range(-2, 3):
                # X=A+sqrt(2)B is monic of degree eight and
                # X(i)=-i(c+sqrt(2)e).
                a = fixed_root_real_part(real([1, -1, 0, 0, 0, 0, 1]), c)
                b = fixed_root_real_part(real([k, 1]), e)
                p0, p1 = add(a, [(0, c)]), add(b, [(0, e)])
                assert at_i(p0) == at_i(p1) == (0, 0)
                norm_p = add(product([p0, p0]), product([p1, p1]), -2)
                expected = scale(add(scale(a, c), scale(b, e), -2), 2)
                assert imag(norm_p) == expected
                assert len(expected) > 1
                count += 1
    # Minimal surviving degree q=1 is possible for a single row.
    # This reaches, rather than exceeds, the d=2 endpoint degree budget.
    a = product([real([1, 0, 1]), real([1, 0, 0, 0, 0, 0, 1])])
    p1 = [(0, -3), (3, 0)]
    assert at_i(a) == at_i(p1) == (0, 0)
    norm_p = add(product([a, a]), product([p1, p1]), -2)
    assert imag(norm_p) == [(0, 0), (36, 0)]
    return count


def check_overlapping_factorization():
    # D includes both t-i and t+i. The residuals also overlap D and
    # each other. Direct product factoring must retain multiplicities.
    minus_i, plus_i = [(0, -1), (1, 0)], [(0, 1), (1, 0)]
    d = product([minus_i, minus_i, plus_i])
    u = product([plus_i, [(2, 3), (1, 0)]])
    v = product([plus_i, [(5, -2), (1, 0)]])
    qi, qj = product([d, u]), product([d, v])
    left = imag(product([conj(qi), qj]))
    right = product([d, conj(d), imag(product([conj(u), v]))])
    assert left == right and any(a for a, _ in left)


def check_quadratic_endpoint_gcd():
    quadratic, t = real([1, 0, 1]), real([0, 1])
    count = 0
    for degree in (4, 6, 8):
        monomial = real([0] * (degree - 2) + [1])
        quotients = [monomial, add(monomial, real([1])),
                     product([quadratic] * ((degree - 2) // 2))]
        for quotient in quotients:
            a = product([quadratic, quotient])
            for d in (2, 3, 5):
                for s in (Fraction(1, 2), 2, -1, -3):
                    f = add(product([a, a]), scale(quadratic, d * s * s))
                    g = scale(product([t, a]), 2 * s)
                    assert f[0][0] > 0
                    assert rational_gcd(f, g) == quadratic
                    count += 1
    return count


def check_cubic_fixture():
    a = product([real([1, 0, 1]), real([4, 0, 1])])
    u, c = real([0, -3, 0, -1]), real([1, 0, 1])
    b = add(u, [(0, 2)])
    assert at_i(a) == at_i(b) == at_i(c) == (0, 0)
    assert product([a, c]) == add(product([u, u]), real([4]))
    # Norm(a+alpha*b+alpha^2*c), alpha^3=2.
    q = add(product([a, a, a]), product([b, b, b]), 2)
    q = add(q, product([c, c, c]), 4)
    q = add(q, product([a, b, c]), -6)
    assert len(q) == 13 and q[-1] == (1, 0)
    assert imag(q) == [(-64, 0)]
    assert -64 == -4 * (2**3 * 2)
    assert at_i(q) == (0, 0)
    assert at_i(derivative(q)) == (0, 0)
    assert at_i(derivative(derivative(q))) == (0, 0)
    assert at_i(derivative(derivative(derivative(q)))) != (0, 0)
    # alpha=3 and i=2 are simple roots modulo 5, so both lift to Q_5.
    assert (3**3 - 2) % 5 == (2**2 + 1) % 5 == 0
    assert (3 * 3**2) % 5 != 0 and (2 * 2) % 5 != 0
    p_mod = add(a, b, 3)
    p_mod = add(p_mod, c, 3**2)
    reduced = [(x + 2 * y) % 5 for x, y in p_mod]
    factored = product([real([0, 1]), real([-2, 1]), real([2, 4, 1])])
    assert reduced == [x % 5 for x, _ in factored]
    assert (4**2 - 4 * 2) % 5 == 3
    assert 3 not in {x * x % 5 for x in range(5)}
    for exponent in range(2, 7):
        # Enlarging the real coefficient field raises this norm to a
        # proper power. Its imaginary polynomial then has positive degree.
        assert len(imag(product([q] * exponent))) > 1


def check_proper_power_degree():
    count = 0
    for n in range(2, 7):
        a = real([7, -3] + [0] * (n - 2) + [1])
        for h in range(n):
            b = [2] if h == 0 else [1] + [0] * (h - 1) + [2]
            q = add(a, [(0, value) for value in b])
            for exponent in range(2, 7):
                imaginary = imag(product([q] * exponent))
                assert len(imaginary) - 1 == (exponent - 1) * n + h
                assert imaginary[-1] == (2 * exponent, 0)
                count += 1
    assert count == 100


def check_constant_brackets():
    blocks = {key: [value, (1, 0)] for key, value in FACTORS.items()}
    rows = [product(blocks[key] for key in incident) for incident in INCIDENT]
    for i, j in combinations(range(3), 2):
        ci, cj = rows[i][0][1], rows[j][0][1]
        assert ci != 0 and cj != 0 and ci != cj
        assert all(b == 0 for _, b in rows[i][1:] + rows[j][1:])
        shared = set(INCIDENT[i]) & set(INCIDENT[j])
        d = product(blocks[key] for key in shared)
        bracket = imag(product([conj(rows[i]), rows[j]]))
        assert bracket == scale(product([d, conj(d)]), cj - ci)
    # Equal nonzero residuals and a half-degree common factor can occur
    # for identical rows; the correct formula allows the zero bracket.
    assert imag(product([conj(rows[0]), rows[0]])) == [(0, 0)]


def main():
    check_parity()
    count = check_quadratic()
    gcd_count = check_quadratic_endpoint_gcd()
    check_overlapping_factorization()
    check_cubic_fixture()
    check_proper_power_degree()
    check_constant_brackets()
    print(f'PASS: 20 parity signs, 11 fixed-root products, {count} quadratic '
          f'norm fixtures, {gcd_count} quadratic endpoint gcds, '
          'repeated/conjugate overlap factoring, the cubic '
          'norm fixture, 100 proper-power degree checks, and three '
          'constant-imaginary pair identities.')


if __name__ == '__main__':
    main()
