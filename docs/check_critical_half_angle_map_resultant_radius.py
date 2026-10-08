"""Exact fixtures for homogeneous half-angle-map content and prime overlap."""

from math import gcd


def mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0]*z[0]+z[1]*z[1]


def nearest(n, d):
    return (n+d//2)//d if n >= 0 else -nearest(-n, d)


def quotient_remainder(z, w):
    denominator = norm(w)
    numerator = mul(z, conj(w))
    q = (nearest(numerator[0], denominator),
         nearest(numerator[1], denominator))
    return q, (z[0]-mul(q, w)[0], z[1]-mul(q, w)[1])


def gaussian_gcd(z, w):
    while w != (0, 0):
        _, remainder = quotient_remainder(z, w)
        z, w = w, remainder
    return z


def gaussian_divides(divisor, value):
    if divisor == (0, 0):
        return value == (0, 0)
    _, remainder = quotient_remainder(value, divisor)
    return remainder == (0, 0)


def pure_power(a, b, degree):
    result = (1, 0)
    for _ in range(degree):
        result = mul(result, (a, b))
    return result


def main():
    pairs = [(a, b) for a in range(1, 120) for b in range(1, 16)
             if gcd(a, b) == 1]
    for a, b in pairs:
        source_denominator = (a, -b)
        F, G = a*a+b*b, b*b
        raw_target_denominator = (F, -G)
        assert gcd(F, G) == 1  # Resultant is one for these two forms.
        assert norm(gaussian_gcd(source_denominator,
                                 raw_target_denominator)) == 1
        # The Gaussian remainder identity is T-b^2(-i)=a^2+b^2=D*bar D.
        assert gaussian_divides(source_denominator, (F, 0))
        assert norm((F, G)) == (a*a+b*b)**2+b**4
        for degree in (2, 3, 5):
            P, Q = pure_power(a, b, degree)
            assert P*P+Q*Q == (a*a+b*b)**degree
            common = gcd(P, Q)
            assert common & (common-1) == 0  # Only the ramified prime can cancel.
        # An order-two cubic can retain the old denominator yet creates
        # a new norm factor. Its conjugate output has D as a divisor.
        cubic_F = a**3+a*b*b-b**3
        cubic_G = a*b*b
        assert gaussian_divides(source_denominator, (cubic_F, -cubic_G))
        assert (cubic_F*cubic_F+cubic_G*cubic_G) % (a*a+b*b) == 0
    print("PASS:", len(pairs), "primitive half-angle pairs: exact degree-two source-prime erasure")
    print("PASS: power-map norm identities, and order-two cubic persistence with a new norm factor")


if __name__ == "__main__":
    main()
