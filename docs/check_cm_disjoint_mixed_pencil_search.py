"""Exact bounded search for a j=1728 disjoint mixed quadratic pencil."""

from fractions import Fraction as F
from itertools import combinations


BOUND = 15


def pencil_data(roots):
    a, b, c, d = map(F, roots)
    r3 = -(a + b + c + d)
    r2 = a*b + a*c + a*d + b*c + b*d + c*d
    r1 = -(a*b*c + a*b*d + a*c*d + b*c*d)
    r0 = a*b*c*d
    p0, p1 = r0, r1 + r0
    q1, q0 = r3, r2 + p1
    if p1 == 0:
        return None
    lambda_b = 4*p0/(p1*p1)
    u, v = p1 - 4*q1, p0 - 4*q0
    if u*u + 16*v == 0:
        return None
    lambda_c = (8*u + 4*v - 64)/(u*u + 16*v)
    if lambda_b == 0 or lambda_c in (0, lambda_b):
        return None

    def discriminant_coefficients(lam):
        d1 = (-2*lam*p1*(1 + lam*q1) + 4*(-1 + lam*q0)
              + 4*lam*lam*p0)
        d2 = (1 + lam*q1)**2 - 4*lam*(-1 + lam*q0)
        return d1, d2

    d1, d2 = discriminant_coefficients(lambda_b)
    if d2 == 0:
        return None
    s = -d1/d2
    d1, d2 = discriminant_coefficients(lambda_c)
    if d2 == 0:
        return None
    t = -d1/d2 - 4
    return (p0, p1, q0, q1, lambda_b, lambda_c, s, t)


def common_value(root):
    return root*root/(root - 1)


def harmonic_equations(s, t):
    return (s*t - 2*s - 2*t,
            s*t + 4*s - 8*t,
            s*t - 8*s + 4*t)


def check_cross_product_identity():
    roots = tuple(map(F, (-3, -2, 4, 5)))
    p0, p1, q0, q1, *_ = pencil_data(roots)
    for x in map(F, range(-8, 9)):
        left = x*x*(x*x + q1*x + q0) - (p1*x + p0)*(x - 1)
        right = F(1)
        for root in roots:
            right *= x - root
        assert left == right


def bounded_search():
    values = [F(n) for n in range(-BOUND, BOUND + 1)
              if n not in (0, 1, 2)]
    admissible = 0
    for roots in combinations(values, 4):
        targets = tuple(common_value(root) for root in roots)
        if len(set(targets)) != 4 or set(targets) & {F(0), F(4)}:
            continue
        data = pencil_data(roots)
        if data is None:
            continue
        *_, s, t = data
        if len({F(0), F(4), s, t}) != 4:
            continue
        admissible += 1
        assert all(value != 0 for value in harmonic_equations(s, t))
    assert admissible == 20452


def check_harmonic_equivalence():
    examples = ((F(1), F(-2), 0),
                (F(-4), F(4, 3), 0))
    for s, t, equation_index in examples:
        cross_ratio = (-s)*(4 - t)/((-t)*(4 - s))
        assert cross_ratio == -1
        assert harmonic_equations(s, t)[equation_index] == 0


def legendre_j(lam):
    return 256*(1 - lam + lam*lam)**3/(lam*lam*(1 - lam)**2)


def cubic_j(a, b):
    """j of y^2=x^3+a*x^2+b*x."""
    return 256*(a*a - 3*b)**3/(b*b*(a*a - 4*b))


def check_common_root_quartic_is_not_a_two_isogeny():
    # The branch values and common arguments of the known disjoint example.
    u, v = F(0), F(4)
    s = F(212248804, 67144321)
    t = F(208624, 46513)
    lam = (u - s)*(v - t)/((u - t)*(v - s))
    assert lam == F(-11090000009, 27202483360)

    # Quotient y^2=x^3+A*x^2+B*x by (0,0):
    # y^2=x^3-2A*x^2+(A^2-4B)*x.  Translate each of the three
    # rational two-torsion points of the Legendre model to zero.
    source_models = ((-(1 + lam), lam),
                     (2 - lam, 1 - lam),
                     (2*lam - 1, lam*(lam - 1)))
    quotient_j = tuple(cubic_j(-2*a, a*a - 4*b)
                       for a, b in source_models)

    roots = tuple(map(F, (-1, 3, 4, 5)))
    lam_r = ((roots[0] - roots[2])*(roots[1] - roots[3])
             / ((roots[0] - roots[3])*(roots[1] - roots[2])))
    assert lam_r == F(5, 3)
    common_root_j = legendre_j(lam_r)
    assert common_root_j == F(438976, 225)
    assert all(value != common_root_j for value in quotient_j)

    # A compact exact certificate for the three inequalities.
    def residue(value, prime):
        return value.numerator*pow(value.denominator, -1, prime) % prime

    assert residue(common_root_j, 101) == 54
    assert tuple(residue(value, 101) for value in quotient_j) == (79, 1, 53)


if __name__ == "__main__":
    check_cross_product_identity()
    check_harmonic_equivalence()
    check_common_root_quartic_is_not_a_two_isogeny()
    bounded_search()
    print("The pencil identity and the three exact harmonic equations pass.")
    print("The known common-root quartic is not a rational 2-isogeny quotient.")
    print("No admissible j=1728 member has four integral common roots in the bounded box.")
