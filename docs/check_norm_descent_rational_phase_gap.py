"""Exact arithmetic checks for norm_descent_rational_phase_gap.md.

The general theorems are proved in the note. These checks independently
exercise the ramified valuations, Gaussian height accounting, and a
nonempty cubic-field example. Standard library only.
"""

from fractions import Fraction as F
from math import gcd


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def neg(z):
    return -z[0], -z[1]


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def scale(z, a):
    return a * z[0], a * z[1]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def power(z, n):
    answer = (1, 0)
    for _ in range(n):
        answer = mul(answer, z)
    return answer


def pi_valuation(z):
    assert z != (0, 0)
    value = 0
    while (z[0] - z[1]) % 2 == 0:
        z = ((z[0] + z[1]) // 2, (z[1] - z[0]) // 2)
        value += 1
    return value


def gaussian_gcd(z, w):
    while w != (0, 0):
        numerator = mul(z, conj(w))
        denominator = norm(w)
        quotient = tuple((2 * a + denominator) // (2 * denominator)
                         for a in numerator)
        remainder = add(z, neg(mul(quotient, w)))
        assert norm(remainder) < norm(w)
        z, w = w, remainder
    return z


def height_squared(p, q):
    assert q != (0, 0)
    divisor_norm = norm(gaussian_gcd(p, q))
    assert norm(p) % divisor_norm == norm(q) % divisor_norm == 0
    return max(norm(p), norm(q)) // divisor_norm


def check_ramified_valuations():
    cases = 0
    for a in range(-12, 13):
        for b in range(1, 13):
            if gcd(a, b) != 1:
                continue
            x, y = (a, b), (a, -b)
            difference = add(x, neg(y))
            assert difference == (0, 2 * b)
            # Degree one is deliberately excluded by the theorem.
            assert difference == power((0, 2 * b), 1)
            for d in range(3, 26, 2):
                lhs = add(power(x, d), neg(power(y, d)))
                rhs = power((0, 2 * b), d)
                if (a - b) % 2:
                    s = pi_valuation(difference)
                    assert pi_valuation(lhs) == s
                    assert pi_valuation(rhs) == d * s
                else:
                    assert a % 2 == b % 2 == 1
                    assert pi_valuation(x) == pi_valuation(y) == 1
                    assert pi_valuation(lhs) == d + 1
                    assert pi_valuation(rhs) == 2 * d
                assert lhs != rhs
                cases += 1
    return cases


def poly_add(a, b):
    return [add(a[j] if j < len(a) else (0, 0),
                b[j] if j < len(b) else (0, 0))
            for j in range(max(len(a), len(b)))]


def poly_mul(a, b):
    result = [(0, 0)] * (len(a) + len(b) - 1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            result[j + k] = add(result[j + k], mul(x, y))
    return result


def poly_scale(a, scalar):
    return [scale(x, scalar) for x in a]


def check_cubic_example():
    # In Q(alpha), alpha^3=2:
    # Norm(x+y*alpha+z*alpha^2)=x^3+2y^3+4z^3-6xyz.
    # The chosen P has x=(t^2+1)(t^2+4),
    # y=-t(t^2+3)+2i, z=t^2+1.
    x = [(4, 0), (0, 0), (5, 0), (0, 0), (1, 0)]
    y = [(0, 2), (-3, 0), (0, 0), (-1, 0)]
    z = [(1, 0), (0, 0), (1, 0)]
    result = [(0, 0)]
    for polynomial, coefficient in ((x, 1), (y, 2), (z, 4)):
        cube = poly_mul(poly_mul(polynomial, polynomial), polynomial)
        result = poly_add(result, poly_scale(cube, coefficient))
    result = poly_add(result, poly_scale(poly_mul(poly_mul(x, y), z), -6))
    assert result[-1] == (1, 0)
    assert result[0] == (68, -64)
    assert all(coefficient[1] == 0 for coefficient in result[1:])
    # P(i)=0 coefficientwise in the field basis.
    for polynomial in (x, y, z):
        at_i = (0, 0)
        for degree, coefficient in enumerate(polynomial):
            at_i = add(at_i, mul(coefficient, power((0, 1), degree)))
        assert at_i == (0, 0)
    # epsilon=1+alpha+alpha^2 has norm one and inverse alpha-1.
    assert 1 + 2 + 4 - 6 == 1
    raw = [-1, 0, 0, 1]  # (alpha-1)(1+alpha+alpha^2)
    assert (raw[0] + 2 * raw[3], raw[1], raw[2]) == (1, 0, 0)


def check_gaussian_heights():
    # For the fixed algebraic target eta=(1+i*sqrt(2))/(1-i*sqrt(2)),
    # the minimal polynomial is 3z^2+2z+3. Pell approximants give
    # eta_n=(v+i*u)/(v-i*u), with u^2-2v^2=1.
    u, v = 3, 2
    gamma_p, gamma_q = (2, 1), (2, -1)
    gamma_height = height_squared(gamma_p, gamma_q)
    assert gamma_height == 5
    cases = 0
    for _ in range(32):
        assert u * u - 2 * v * v == 1
        p, q = (v, u), (v, -u)
        h_eta = height_squared(p, q)
        assert h_eta == u * u + v * v
        cleared_f = add(add(scale(mul(p, p), 3), scale(mul(p, q), 2)),
                        scale(mul(q, q), 3))
        assert cleared_f == (-4, 0)
        # For A=2 and L=8, the Liouville bound gives H(eta_n)>=v/4.
        assert 16 * h_eta >= v * v
        correction_p, correction_q = mul(gamma_p, q), mul(gamma_q, p)
        h_correction = height_squared(correction_p, correction_q)
        product_p, product_q = mul(correction_p, p), mul(correction_q, q)
        assert mul(product_p, gamma_q) == mul(product_q, gamma_p)
        assert height_squared(product_p, product_q) == gamma_height
        assert h_eta <= h_correction * gamma_height
        assert h_correction <= gamma_height * h_eta
        # The note's general bound with D=2, delta=2, C=4 yields
        # H(correction)^10 >= v^2/80^4. No floating point is used.
        assert 80**4 * h_correction**5 >= v * v
        cases += 1
        u, v = 3 * u + 4 * v, 2 * u + 3 * v
    return cases


def check_algebraic_offset_escape():
    """Exact degree-16 counterexample to removing the rational-offset assumption."""
    one = [(F(1), F(0))]
    g = [(F(x, 5), F(0)) for x in (6, -7, 10, -2, 4)]
    h = [(F(x, 100), F(0)) for x in (-14, -37, 20, -12, 34)]
    b = poly_scale(poly_add(poly_mul(g, g), one), 2)
    a = poly_add(g, poly_mul(b, h))
    c = poly_add(poly_add(one, poly_scale(poly_mul(g, h), 4)),
                 poly_scale(poly_mul(b, poly_mul(h, h)), 2))
    identity = poly_add(poly_add(poly_mul(a, a), one),
                        poly_scale(poly_mul(b, c), F(-1, 2)))
    assert all(z == (0, 0) for z in identity)
    assert (len(a)-1, len(b)-1, len(c)-1) == (12, 8, 16)
    assert c[-1] == (F(4624, 15625), 0)

    def at_i(poly):
        result = (0, 0)
        for j, coefficient in enumerate(poly):
            result = add(result, mul(coefficient, power((0, 1), j)))
        return result

    assert at_i(g) == (0, -1) and at_i(h) == (0, F(-1, 4))
    assert at_i(a) == (0, -1) and at_i(b) == at_i(c) == (0, 0)

    def kmul(x, y):
        raw = [F(0)] * 5
        for j in range(3):
            for k in range(3):
                raw[j+k] += x[j]*y[k]
        return (raw[0]+raw[3]/2, raw[1]+raw[4]/2, raw[2])

    def kadd(x, y):
        return tuple(a+b for a, b in zip(x, y))

    def at_alpha(poly):
        result = (F(0), F(0), F(0))
        alpha = (F(0), F(1), F(0))
        for coefficient in reversed(poly):
            assert coefficient[1] == 0
            result = kadd(kmul(result, alpha), (coefficient[0], F(0), F(0)))
        return result

    assert at_alpha(g) == (1, -1, 2)
    assert at_alpha(h) == (F(-1, 5), F(-1, 5), F(1, 5))
    assert at_alpha(a) == (0, 0, 0)
    assert at_alpha(b) == (0, 0, 10)
    assert at_alpha(c) == (0, F(2, 5), 0)
    alpha = (F(0), F(1), F(0))
    alpha_squared = (F(0), F(0), F(1))
    u_value = kadd(at_alpha(a),
                   kadd(kmul(at_alpha(b), alpha), kmul(at_alpha(c), alpha_squared)))
    assert u_value == (F(26, 5), 0, 0)

    # Contact order in cubic_pisot_norm_path_endpoint_obstruction.md is one.
    def derivative(poly):
        return [scale(poly[j], j) for j in range(1, len(poly))]

    u_derivative = kadd(at_alpha(derivative(a)), kadd(
        kmul(at_alpha(derivative(b)), alpha),
        kmul(at_alpha(derivative(c)), alpha_squared)))
    assert u_derivative == (F(-987, 125), F(1599, 125), F(4371, 250))
    assert u_derivative != (0, 0, 0)

    # Norm(U+i) for alpha^3=1/2; imaginary part is exactly -4.
    x = poly_add(a, [(0, 1)])
    norm_poly = [(0, 0)]
    for polynomial, coefficient in ((x, 1), (b, F(1, 2)), (c, F(1, 4))):
        norm_poly = poly_add(norm_poly, poly_scale(
            poly_mul(poly_mul(polynomial, polynomial), polynomial), coefficient))
    norm_poly = poly_add(norm_poly, poly_scale(poly_mul(poly_mul(x, b), c), F(-3, 2)))
    assert norm_poly[0][1] == -4
    assert all(z[1] == 0 for z in norm_poly[1:])
    leading_u = (F(0), F(0), F(4624, 15625))
    y = (F(0), F(15625, 2312), F(0))
    assert kmul(leading_u, y) == (1, 0, 0)
    epsilon = (F(1), F(2), F(2))
    epsilon_inverse = (F(-1), F(0), F(2))
    assert kmul(epsilon, epsilon_inverse) == (1, 0, 0)
    beta = (F(0), F(0), F(2))
    assert kmul(kmul(beta, beta), beta) == (2, 0, 0)
    assert kadd(kadd((F(1), F(0), F(0)), beta), kmul(beta, beta)) == epsilon


if __name__ == "__main__":
    ramified = check_ramified_valuations()
    check_cubic_example()
    heights = check_gaussian_heights()
    check_algebraic_offset_escape()
    print(f"PASS: {ramified} odd-power ramified valuation cases;")
    print("      exact cubic polynomial norm, fixed root, and unit identities;")
    print(f"      {heights} Pell phase and Gaussian-height fixtures.")
    print("PASS: degree-16 offset, contact order one, exact norm and unit identities.")
