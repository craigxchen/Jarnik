"""Exact certificates for the rational all-split q=2 norm-contact frame.

Only the Python standard library is needed. Polynomial coefficients are
stored in increasing powers of t. Finite-field certificates establish
squarefreeness, coprimality, and irreducibility over the rationals.
"""

from itertools import combinations


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    return trim([
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ])


def scale(a, c):
    return trim([c * x for x in a])


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))] or [0])


def gcd_mod(a, b, prime=7):
    a = trim([x % prime for x in a])
    b = trim([x % prime for x in b])
    while b != [0]:
        r = a[:]
        while r != [0] and len(r) >= len(b):
            c = r[-1] * pow(b[-1], -1, prime) % prime
            shift = len(r) - len(b)
            for i, x in enumerate(b):
                r[shift + i] = (r[shift + i] - c * x) % prime
            r = trim(r)
        a, b = b, r
    inv = pow(a[-1], -1, prime)
    return [(x * inv) % prime for x in a]


def norm_numerator(real, imag):
    return add(mul(real, real), mul(imag, imag))


def remainder_mod(a, modulus, prime):
    a = trim([x % prime for x in a])
    modulus = trim([x % prime for x in modulus])
    inverse_lead = pow(modulus[-1], -1, prime)
    while a != [0] and len(a) >= len(modulus):
        c = a[-1] * inverse_lead % prime
        shift = len(a) - len(modulus)
        for i, x in enumerate(modulus):
            a[shift + i] = (a[shift + i] - c * x) % prime
        a = trim(a)
    return a


def power_mod(a, exponent, modulus, prime):
    out = [1]
    a = remainder_mod(a, modulus, prime)
    while exponent:
        if exponent & 1:
            out = remainder_mod(mul(out, a), modulus, prime)
        a = remainder_mod(mul(a, a), modulus, prime)
        exponent //= 2
    return out


def certify_irreducible_degree_14(n, prime):
    """Rabin's criterion; the prime divisors of degree 14 are 2 and 7."""
    assert prime >= 2
    assert all(prime % d for d in range(2, int(prime**0.5) + 1))
    assert len(n) - 1 == 14 and n[-1] % prime != 0
    x = [0, 1]
    frobenius = x
    for k in range(1, 15):
        frobenius = power_mod(frobenius, prime, n, prime)
        if k in (2, 7):
            assert gcd_mod(n, add(frobenius, scale(x, -1)), prime) == [1]
    assert frobenius == x


def main():
    q = [1, 0, 0, 0, 0, 0, 0, 1]
    s = [0, 0, 0, -1, 1]
    p = add(s, scale(q, -2))
    p_minus_q = add(p, scale(q, -1))
    two_p_minus_q = add(scale(p, 2), scale(q, -1))
    p_plus_two_q = add(p, scale(q, 2))
    numerators = [
        norm_numerator(p, q),
        norm_numerator(p_minus_q, q),
        norm_numerator(two_p_minus_q, scale(q, 2)),
        norm_numerator(p_plus_two_q, q),
    ]
    g = [1] + [0] * 7 + [1]
    h = [1] + [0] * 5 + [1]

    # The frame determinant is (P-(P-Q))/Q = 1.
    assert add(p, scale(p_minus_q, -1)) == q
    assert p_plus_two_q == s
    assert numerators[3] == mul(g, h)
    assert len(p) - 1 == len(q) - 1 == 7
    assert gcd_mod(p, q) == [1]
    assert gcd_mod(q, derivative(q), 11) == [1]
    assert [len(n) - 1 for n in numerators] == [14] * 4
    assert [n[-1] for n in numerators] == [5, 10, 29, 1]
    for n in numerators:
        assert n[-1] % 7 != 0
        assert gcd_mod(n, q) == [1]
        assert gcd_mod(n, derivative(n)) == [1]

    # Pairwise disjointness is proved by the distinct target values in
    # the note. Use exact mod 11 certificates as an independent check;
    # the target values can merge after reduction mod 7.
    for a, b in combinations(numerators, 2):
        assert a[-1] % 11 and b[-1] % 11
        assert gcd_mod(a, b, 11) == [1]
    assert gcd_mod(g, q) == [1]
    assert gcd_mod(g, derivative(g)) == [1]
    assert [len(n) - 1 for n in numerators[:3]] + [len(g) - 1] == [14, 14, 14, 8]
    for n, prime in zip(numerators[:3], (43, 127, 59)):
        certify_irreducible_degree_14(n, prime)
    print("PASS: determinant one, pole degree 7, contact degrees 14,14,14,8.")
    print("PASS: exact factorization; all contacts reduced, pole-free, and disjoint.")
    print("PASS: N1,N2,N3 irreducible modulo 43,127,59 respectively, hence over Q.")
    print("Thus their degree-14 divisors cannot split into Q-defined degree-2 contacts.")
    print("Split residue fields follow from the norm identities, as proved in the note.")


if __name__ == "__main__":
    main()
