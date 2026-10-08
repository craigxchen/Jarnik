"""Exact finite inputs for the all-six-place 17-adic escape proposition."""

from fractions import Fraction as Q


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, c):
    return tuple(c*x for x in a)


def mul(a, b):
    raw = [Q(0)] * 5
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            raw[j+k] += x*y
    # alpha^3=3alpha-1 and alpha^4=3alpha^2-alpha.
    return (raw[0]-raw[3], raw[1]+3*raw[3]-raw[4], raw[2]+3*raw[4])


def reduce_fraction(x, modulus):
    x = Q(x)
    return x.numerator * pow(x.denominator, -1, modulus) % modulus


def evaluate(a, alpha, modulus):
    return sum(reduce_fraction(x, modulus) * pow(alpha, j, modulus)
               for j, x in enumerate(a)) % modulus


def lift_cubic_root(a):
    assert (a*a*a-3*a+1) % 17 == 0
    derivative = (3*a*a-3) % 17
    assert derivative
    correction = -((a*a*a-3*a+1)//17) * pow(derivative, -1, 17) % 17
    lifted = a+17*correction
    assert (lifted**3-3*lifted+1) % (17*17) == 0
    return lifted


def hensel_lift(a, prime, exponent, polynomial, derivative):
    modulus = prime
    assert polynomial(a) % modulus == 0
    for _ in range(1, exponent):
        correction = -(polynomial(a) // modulus) * pow(derivative(a), -1, prime)
        a += modulus * (correction % prime)
        modulus *= prime
        assert polynomial(a) % modulus == 0
    return a


def residual_coefficients(alpha, imaginary_root, rotation_p, rotation_q, r):
    """Ascending monic cubic coefficients from the analytic formula (1)."""
    w, v = 2-alpha*alpha, alpha*alpha+2*alpha-2
    W = rotation_p*w+rotation_q*v
    V = -3*rotation_q*w+rotation_p*v
    L = 1+W
    u = 2*r/(1+3*r*r)
    s = (1-3*r*r)/(1+3*r*r)
    h = 3*u*u
    A, B = 2*h/s, 2*u/s
    i = imaginary_root
    b = i*L-A*L+2*u*V
    c = 1+W*(1-2*h)-i*L*A+V*B*(-h-1+i*s)
    d = -2*A+i*(1+W*(1-2*h)-V*B*(h+1))
    return (d/L, c/L, b/L, Q(1))


def discriminant(coefficients):
    d, c, b, a = coefficients
    return b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d


def rational_valuation_unit(x, prime):
    assert x
    numerator, denominator = x.numerator, x.denominator
    valuation = 0
    while numerator % prime == 0:
        numerator //= prime
        valuation += 1
    while denominator % prime == 0:
        denominator //= prime
        valuation -= 1
    return valuation, numerator * pow(denominator, -1, prime) % prime


def check_parameter_fixtures():
    # Comparing a rational coefficient at two lifts can lose two powers
    # of 17 through its pair of L denominators. Twelve-digit lifts still
    # determine each coefficient mod 17^10. Even allowing another three
    # denominator powers in the degree-four discriminant gives precision
    # mod 17^7, enough for valuations <=6 and their unit parts.
    for i_residue in (4, 13):
        ii = hensel_lift(i_residue, 17, 12, lambda x: x*x+1, lambda x: 2*x)
        results = []
        for a_residue in (7, 13, 14):
            aa = hensel_lift(a_residue, 17, 12,
                             lambda x: x**3-3*x+1, lambda x: 3*x*x-3)
            coef = residual_coefficients(Q(aa), Q(ii), Q(-13, 14), Q(-3, 14), Q(17**3))
            assert all(rational_valuation_unit(x, 17)[0] >= 0 for x in coef if x)
            results.append(rational_valuation_unit(discriminant(coef), 17))
        expected_unit = 8 if i_residue == 4 else 9
        assert results == [(4, expected_unit), (6, expected_unit), (4, expected_unit)]

    # A separate prime in the same class 17 mod 36 has six simple split
    # reductions. This is a local fixture, not global splitting over K(i).
    prime = 197
    p, q, r = Q(-971, 973), Q(36, 973), Q(102)
    assert p*p+3*q*q == 1 and prime % 36 == 17
    for aa in (34, 169, 191):
        assert (aa**3-3*aa+1) % prime == 0
        for ii in (14, 183):
            assert (ii*ii+1) % prime == 0
            coef = tuple(reduce_fraction(x, prime) for x in
                         residual_coefficients(Q(aa), Q(ii), p, q, r))
            roots = [t for t in range(prime)
                     if sum(c*pow(t, j, prime) for j, c in enumerate(coef)) % prime == 0]
            assert len(roots) == 3 and discriminant(coef) % prime


def main():
    one = (Q(1), Q(0), Q(0))
    w = (Q(2), Q(0), Q(-1))
    v = (Q(-2), Q(2), Q(1))
    p, q = Q(-13, 14), Q(-3, 14)
    assert p*p+3*q*q == 1
    W = add(scale(w, p), scale(v, q))
    V = add(scale(w, -3*q), scale(v, p))
    L = add(one, W)
    assert W == (Q(-10, 7), Q(-3, 7), Q(5, 7))
    assert V == (Q(22, 7), Q(-13, 7), Q(-11, 7))
    assert L == (Q(-3, 7), Q(-3, 7), Q(5, 7))
    assert add(mul(W, W), scale(mul(V, V), Q(1, 3))) == scale(one, 4)

    roots = (7, 13, 14)
    lifts = tuple(map(lift_cubic_root, roots))
    assert lifts == (75, 132, 82)
    leading_residues = tuple(evaluate(L, a, 17) for a in roots)
    assert leading_residues == (0, 3, 0)
    e_values = []
    for a in lifts:
        lead = evaluate(L, a, 17*17)
        assert lead  # Nonzero mod 17^2 proves e <= 1.
        e_values.append(int(lead % 17 == 0))
    assert e_values == [1, 0, 1]
    assert evaluate(L, lifts[0], 17*17)//17 == 11
    assert evaluate(L, lifts[2], 17*17)//17 == 6

    for imaginary_root in (4, 13):
        assert (imaginary_root**2+1) % 17 == 0
        assert 2*imaginary_root % 17
        numerators = []
        for a in roots:
            lead = evaluate(L, a, 17)
            vv = evaluate(V, a, 17)
            numerator = -128*imaginary_root*(vv*vv-12*lead) % 17
            assert numerator and pow(numerator, 8, 17) == 1
            numerators.append(numerator)
        assert numerators == ([16, 4, 16] if imaginary_root == 4 else [1, 13, 1])

    # The proof uses k-4e > -2e for every k>=3 and e in {0,1}.
    assert all(3-4*e > -2*e for e in e_values)
    check_parameter_fixtures()
    print('PASS: exact rational rotation and cubic-field identity.')
    print('PASS: three Hensel lifts, lead valuations (1,0,1), and all six square residues.')
    print('PASS: six discriminant fixtures at r=17^3 and six simple split reductions at 197.')
    print('The analytic remainder argument proves local splitting for every v_17(r)>=3.')


if __name__ == '__main__':
    main()
