"""Exact checks for disjoint_mixed_elliptic_six_curve.md."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt


LAMBDAS = (F(0), F(-240, 5329), F(8, 257))
COMMON_X = (F(-1), F(3), F(4), F(5))
COMMON_TARGETS = (F(-1, 2), F(9, 2), F(16, 3), F(25, 4))
VARIABLES = ("a", "b", "c")
LABELS = frozenset(("inf", "0", "1") + VARIABLES)


def coefficients(lam):
    # Low-to-high coefficients.
    numerator = (-60 * lam, -73 * lam, F(1))
    denominator = (-1 - 38 * lam, 1 - 11 * lam, lam)
    return numerator, denominator


def value(lam, x):
    numerator, denominator = coefficients(lam)
    top = sum(coef * x**index for index, coef in enumerate(numerator))
    bottom = sum(coef * x**index for index, coef in enumerate(denominator))
    assert bottom
    return top / bottom


def critical_data(lam):
    numerator, denominator = coefficients(lam)
    p0, p1, p2 = numerator
    q0, q1, q2 = denominator
    derivative = (
        p1 * q0 - p0 * q1,
        2 * (p2 * q0 - p0 * q2),
        p2 * q1 - p1 * q2,
    )
    c0, c1, c2 = derivative
    discriminant = c1 * c1 - 4 * c2 * c0
    root_num = isqrt(discriminant.numerator)
    root_den = isqrt(discriminant.denominator)
    assert root_num**2 == discriminant.numerator
    assert root_den**2 == discriminant.denominator
    square_root = F(root_num, root_den)
    roots = ((-c1 + square_root) / (2 * c2),
             (-c1 - square_root) / (2 * c2))
    return roots, tuple(value(lam, root) for root in roots)


def canonical_boundary(side):
    side = frozenset(side)
    other = LABELS - side
    assert len(side) >= 2 and len(other) >= 2
    if len(side) != len(other):
        return side if len(side) < len(other) else other
    return min(side, other, key=lambda part: tuple(sorted(part)))


def all_boundaries():
    result = set()
    for size in (2, 3, 4):
        for side in combinations(sorted(LABELS), size):
            if len(LABELS - set(side)) >= 2:
                result.add(canonical_boundary(side))
    return result


def check_pencil():
    for x, target in zip(COMMON_X, COMMON_TARGETS):
        assert tuple(value(lam, x) for lam in LAMBDAS) == (target,) * 3

    expected = (
        ((F(2), F(0)), (F(4), F(0))),
        ((F(402794, 639337), F(-120, 73)),
         (F(212248804, 67144321), F(0))),
        ((F(3244, 1069), F(14, 5)),
         (F(208624, 46513), F(4))),
    )
    actual = tuple(critical_data(lam) for lam in LAMBDAS)
    assert actual == expected
    branch_sets = tuple(set(values) for _, values in actual)
    assert branch_sets[0] & branch_sets[1] == {F(0)}
    assert branch_sets[0] & branch_sets[2] == {F(4)}
    assert not branch_sets[1] & branch_sets[2]
    assert len(set().union(*branch_sets)) == 4
    assert not set(COMMON_TARGETS) & set().union(*branch_sets)

    # Exact affine determinant; its four simple roots are COMMON_X.
    for lam, mu in combinations(LAMBDAS, 2):
        for x in map(F, range(-6, 8)):
            p_l, q_l = coefficients(lam)
            p_m, q_m = coefficients(mu)
            pl = sum(coef * x**i for i, coef in enumerate(p_l))
            ql = sum(coef * x**i for i, coef in enumerate(q_l))
            pm = sum(coef * x**i for i, coef in enumerate(p_m))
            qm = sum(coef * x**i for i, coef in enumerate(q_m))
            determinant = pl * qm - pm * ql
            expected_det = ((mu - lam) * (x + 1) * (x - 3)
                            * (x - 4) * (x - 5))
            assert determinant == expected_det

    # The common source normalization.
    transform = lambda X: (16 * X - 1) / (4 * X + 1)
    assert tuple(transform(X) for X in (F(0), F(1), F(-3, 2))) == (
        F(-1), F(3), F(5))
    # At X=infinity the ratio of leading coefficients is 4.
    assert F(16, 4) == 4


def check_disjoint_boundaries():
    contacts = Counter()
    for anchor in ("0", "1", "inf"):
        for size in (1, 2, 3):
            for subset in combinations(VARIABLES, size):
                contacts[(anchor, subset)] += 1
    for size in (2, 3):
        for subset in combinations(VARIABLES, size):
            contacts[("rho", subset)] += 1
    assert len(contacts) == 25
    assert set(contacts.values()) == {1}

    divisors = {
        canonical_boundary({anchor} | set(subset))
        for anchor, subset in contacts if anchor != "rho"
    } | {
        canonical_boundary(subset)
        for anchor, subset in contacts if anchor == "rho"
    }
    assert divisors == all_boundaries()
    assert len(divisors) == 25


def _poly_add(*polynomials):
    length = max(map(len, polynomials))
    result = [F(0)] * length
    for polynomial in polynomials:
        for index, coefficient in enumerate(polynomial):
            result[index] += coefficient
    return result


def _poly_mul(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    return result


def _transform_coefficients(coefficients):
    # If x=(16X-1)/(4X+1), clear the common denominator in
    # p(x) by replacing (1,x,x^2) with (h^2,nh,n^2).
    h = (F(1), F(4))
    n = (F(-1), F(16))
    return _poly_add(
        [coefficients[0] * coefficient for coefficient in _poly_mul(h, h)],
        [coefficients[1] * coefficient for coefficient in _poly_mul(n, h)],
        [coefficients[2] * coefficient for coefficient in _poly_mul(n, n)],
    )


def _linear_t_mul(left, right):
    # Multiplication of two degree-at-most-one polynomials in T.
    return (
        left[0] * right[0],
        left[0] * right[1] + left[1] * right[0],
        left[1] * right[1],
    )


def _discriminant_in_T(lam):
    numerator = _transform_coefficients(coefficients(lam)[0])
    denominator = _transform_coefficients(coefficients(lam)[1])
    # The equation in X is N(X)-T D(X)=0.
    a = (numerator[2], -denominator[2])
    b = (numerator[1], -denominator[1])
    c = (numerator[0], -denominator[0])
    bb = _linear_t_mul(b, b)
    ac = _linear_t_mul(a, c)
    return [bb[i] - 4 * ac[i] for i in range(3)]


def _is_rational_square(value):
    if value < 0:
        return False
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    return (numerator_root * numerator_root == value.numerator
            and denominator_root * denominator_root == value.denominator)


def _evaluate_polynomial(coefficients, value):
    return sum(coefficient * value**index
               for index, coefficient in enumerate(coefficients))


def check_quotient_rank_certificate():
    # The three exact discriminants of the inverse quadratics after the
    # common-source normalization.  Keep their original square factors:
    # this is what records the quadratic twist in their product.
    discriminants = tuple(_discriminant_in_T(lam) for lam in LAMBDAS)
    assert discriminants == (
        [F(0), F(-1600), F(400)],
        [F(0), F(-84899521600, 28398241),
         F(26857728400, 28398241)],
        [F(333798400, 66049), F(-157870400, 66049),
         F(18605200, 66049)],
    )

    # Four common unramified fibers split completely over Q.  Each of the
    # three inverse quadratics has two rational roots at each target.
    for target in COMMON_TARGETS:
        values = [_evaluate_polynomial(discriminant, target)
                  for discriminant in discriminants]
        assert all(value and _is_rational_square(value) for value in values)

    # Extract only a rational square from Delta_b Delta_c.  The resulting
    # quartic is the actual twisted model, rather than a monic rescaling.
    qb = [F(0), F(-212248804), F(67144321)]
    qc = [F(834496), F(-394676), F(46513)]
    assert discriminants[1] == [F(400, 73**4) * value for value in qb]
    assert discriminants[2] == [F(400, 257**2) * value for value in qc]
    assert qc == _poly_mul([F(-4), F(1)], [F(-208624), F(46513)])
    quartic = _poly_mul(qb, qc)
    assert quartic[4] == F(67144321 * 46513)
    assert not _is_rational_square(quartic[4])
    witness_t = F(-1, 2)
    witness_y = F(45299439, 4)
    assert _evaluate_polynomial(quartic, witness_t) == witness_y**2
    assert witness_y
    assert witness_t not in (
        F(0), F(212248804, 67144321), F(4), F(208624, 46513))

    def residue(numerator, denominator, prime):
        assert denominator % prime
        return numerator * pow(denominator, -1, prime) % prime

    def point_count(prime):
        assert prime % 2
        roots = [residue(0, 1, prime),
                 residue(212248804, 67144321, prime),
                 residue(4, 1, prime),
                 residue(208624, 46513, prime)]
        assert len(set(roots)) == 4
        squares = {value * value % prime for value in range(prime)}
        finite = 0
        for t in range(prime):
            # Evaluate the integer quartic itself.  Replacing this by the
            # monic root product would silently discard the twist.
            rhs = (t * (67144321 * t - 212248804) * (t - 4)
                   * (46513 * t - 208624)) % prime
            finite += (1 if rhs == 0 else 2 if rhs in squares else 0)
        leading = (67144321 * 46513) % prime
        infinity = 2 if leading in squares else 0
        return finite + infinity, roots, finite, infinity, leading in squares

    first = point_count(23)
    second = point_count(37)
    assert first == (28, [0, 16, 4, 2], 26, 2, True)
    assert second == (40, [0, 11, 4, 23], 40, 0, False)
    torsion_bound = gcd(first[0], second[0])
    assert torsion_bound == 4
    assert torsion_bound < 32


if __name__ == "__main__":
    check_pencil()
    check_disjoint_boundaries()
    check_quotient_rank_certificate()
    print("The rational pencil has the mixed four-branch path exactly.")
    print("Its four common arguments have four distinct target values.")
    print("All 25 stable boundaries occur at distinct source points.")
    print("The twisted E/E[2] quartic has 28 and 40 points mod 23 and 37.")
    print("Its rational torsion order divides gcd(28,40)=4 < 32.")
    print("The quartic point (-1/2,45299439/4) is an exact nontorsion witness.")
