"""Exact certificates for rational specializations of the degree-seven frame."""

from fractions import Fraction
from math import gcd

from check_scaled_split_contact_profile_q2 import (
    add, derivative, gcd_mod, mul, scale,
)


def homogeneous(a, b):
    q = a**7 + b**7
    ss = a**3 * b**3 * (a - b)
    return ss - 2 * q, q


def main():
    count = 0
    integral_parameters = set()
    for a in range(-80, 81):
        for b in range(0, 81):
            if gcd(a, b) != 1 or a + b == 0:
                continue
            p, q = homogeneous(a, b)
            common = gcd(p, q)
            assert common in (1, 2)
            den = abs(q) // common
            hh = max(abs(a), abs(b))
            assert 4 * den >= hh**6
            if den == 1:
                parameter = 'infinity' if b == 0 else Fraction(a, b)
                integral_parameters.add(parameter)
                assert p // q == -2
            count += 1
    assert integral_parameters == {Fraction(0), Fraction(1), 'infinity'}

    for n in range(1, 101):
        p, q = homogeneous(1 - n, n)
        assert gcd(p, q) == 1
        assert q == n**7 - (n - 1)**7

    qpoly = [1] + [0] * 6 + [1]
    spoly = [0, 0, 0, -1, 1]
    d8 = [-3, 4, 0, 0, 0, 0, 0, 4, -3]
    numerator = add(mul(derivative(spoly), qpoly),
                    scale(mul(spoly, derivative(qpoly)), -1))
    assert numerator == [0, 0] + d8
    assert gcd_mod(d8, derivative(d8), 5) == [1]
    assert gcd_mod(d8, qpoly, 5) == [1]
    assert d8[0] != 0
    print(f'PASS: denominator bound and integral-value list on {count} primitive pairs.')
    print('PASS: sharp exponent-six sequence and derivative/fiber certificates.')


if __name__ == '__main__':
    main()
