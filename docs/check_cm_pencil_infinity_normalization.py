"""Exact checks for the common-root-at-infinity CM-pencil normalization.

The search is deliberately a necessary-filter search.  A tuple surviving
the displayed rational denominators is not asserted to give a valid pencil;
the optional degree/gcd and distinct-branch checks are kept separate.
"""

from fractions import Fraction as Q
from itertools import combinations


def f0(x):
    return x * x / (x - 1)


def phi(x, a):
    return a * x / (1 + (a - 1) * x / 2)


def mobius(y, a):
    return a * a * y / (1 + (a * a - 1) * y / 4)


def mobius_inverse(y, a):
    return 4 * y / (4 * a * a - (a * a - 1) * y)


def check_domain_target_normalizer():
    # f0(phi_a(x)) = M_a(f0(x)); all quantities are rational.
    for a in (Q(2), Q(3, 2), Q(-1, 2), Q(5, 3)):
        for x in (Q(-3), Q(-1), Q(1, 2), Q(3), Q(7)):
            if x == 1 or 1 + (a - 1) * x / 2 == 0:
                continue
            assert f0(phi(x, a)) == mobius(f0(x), a)

    # The transformed pencil is f_new=M_a^{-1} o f o phi_a, so f0 stays f0.
    for a in (Q(2), Q(3, 2), Q(-1, 2)):
        for x in (Q(-3), Q(-1), Q(1, 2), Q(3), Q(7)):
            if x == 1 or 1 + (a - 1) * x / 2 == 0:
                continue
            assert mobius_inverse(f0(phi(x, a)), a) == f0(x)

    # If rho is a finite rational common root away from 0,2, then
    # a=rho/(rho-2) gives phi_a(infinity)=rho.  Thus rho is new infinity.
    for rho in (Q(-3), Q(-1), Q(1, 2), Q(3), Q(7)):
        a = rho / (rho - 2)
        assert a != 0
        assert 2 * a / (a - 1) == rho
        # The denominator of phi_a(x)=rho vanishes at the new point infinity.
        assert a - rho * (a - 1) / 2 == 0


def pencil_data(roots):
    r1, r2, r3 = roots
    s1 = r1 + r2 + r3
    s2 = r1 * r2 + r1 * r3 + r2 * r3
    c = -r1 * r2 * r3
    b = s2 + c
    e = b - s1

    # Necessary algebraic denominator filters only.
    if b == 0 or b == 4:
        return None
    lb = 4 * c / (b * b)
    h = 2 * b - 8 + c - 4 * e
    lc = 4 * h / ((b - 4) * (b - 4))
    if c == 0 or h == 0 or lb == lc:
        return None
    if lb == -1 or lc == -1:
        return None

    m = b * b + 4 * c
    k = b * b + 4 * c - 16 * (e + 1)
    if m == 0 or k == 0:
        return None

    s = 4 * b * (b**3 + 2 * c * b * (b - 2 * e) + 8 * c * c) / (m * m)
    t = (8 * h * (b * b * (b - 4 - 2 * e) + 4 * c * (b - 2))) / (k * k)
    harmonic = (
        s * t - 2 * s - 2 * t,
        s * t + 4 * s - 8 * t,
        s * t - 8 * s + 4 * t,
    )
    return b, c, e, lb, lc, s, t, harmonic


def degree_two_map(k, b, c, e):
    """Check that (x^2+k(bx+c))/(x-1+k(x+e)) is degree two and reduced."""
    q1 = 1 + k
    q0 = -1 + k * e
    if q1 == 0 and q0 == 0:
        return False
    # Resultant of a monic quadratic and q1*x+q0.
    resultant = q0 * q0 - k * b * q1 * q0 + k * c * q1 * q1
    return resultant != 0


def search(values):
    candidate_count = 0
    harmonic_hits = []
    admissible_hits = []
    for roots in combinations(values, 3):
        if any(root in (0, 2) for root in roots):
            continue
        data = pencil_data(roots)
        if data is None:
            continue
        candidate_count += 1
        b, c, e, lb, lc, s, t, harmonic = data
        if not any(value == 0 for value in harmonic):
            continue
        harmonic_hits.append((roots, data))
        if (len({Q(0), Q(4), s, t}) == 4
                and degree_two_map(lb, b, c, e)
                and degree_two_map(lc, b, c, e)):
            admissible_hits.append((roots, data))
    return candidate_count, harmonic_hits, admissible_hits


def main():
    check_domain_target_normalizer()

    integers = [Q(x) for x in range(-15, 16) if x not in (0, 2)]
    count, harmonic, admissible = search(integers)
    assert count == 3617
    assert harmonic == []
    assert admissible == []
    print(f"integer candidates={count}; harmonic hits={len(harmonic)}; "
          f"admissible hits={len(admissible)}")

    rationals = sorted({Q(p, q) for q in range(1, 6)
                        for p in range(-15, 16) if Q(p, q) not in (0, 2)})
    count, harmonic, admissible = search(rationals)
    assert count == 186481
    assert harmonic == []
    assert admissible == []
    print(f"rational candidates={count}; harmonic hits={len(harmonic)}; "
          f"admissible hits={len(admissible)}")


if __name__ == "__main__":
    main()
