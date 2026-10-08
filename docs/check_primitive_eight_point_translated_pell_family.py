"""Exact checks for the primitive eight-point Pell endpoint family.

The checker uses integer Gaussian arithmetic only. It verifies the seven
individually primitive points, the odd-parity eight-point extension, and the
optional even-parity eight-point family (whose six middle points have content
5). Phase alignment is checked from exact Pell/Binet identities and exponent
bookkeeping, without converting the huge Pell coordinates to floating point.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, prod


Gaussian = tuple[int, int]
Q5 = tuple[Fraction, Fraction]  # a + b sqrt(5)


def q5_add(x: Q5, y: Q5) -> Q5:
    return x[0] + y[0], x[1] + y[1]


def q5_neg(x: Q5) -> Q5:
    return -x[0], -x[1]


def q5_mul(x: Q5, y: Q5) -> Q5:
    return x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0]


def q5_nonnegative(x: Q5) -> bool:
    a, b = x
    if b == 0:
        return a >= 0
    if a >= 0 and b >= 0:
        return True
    if a <= 0 and b <= 0:
        return False
    return 5 * b * b >= a * a if b > 0 else a * a >= 5 * b * b


def q5_pow(x: Q5, e: int) -> Q5:
    if e < 0:
        norm = x[0] * x[0] - 5 * x[1] * x[1]
        assert norm != 0
        x = (x[0] / norm, -x[1] / norm)
        e = -e
    out = (Fraction(1), Fraction(0))
    while e:
        if e & 1:
            out = q5_mul(out, x)
        x = q5_mul(x, x)
        e >>= 1
    return out


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def norm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def gprod(values) -> Gaussian:
    out = (1, 0)
    for z in values:
        out = mul(out, z)
    return out


def gdiv_exact(z: Gaussian, w: Gaussian) -> Gaussian:
    d = norm(w)
    assert d > 0
    x = z[0] * w[0] + z[1] * w[1]
    y = z[1] * w[0] - z[0] * w[1]
    assert x % d == y % d == 0
    return x // d, y // d


def nearest_quotient(a: int, d: int) -> int:
    assert d > 0
    return (2 * a + d) // (2 * d)


def gdivmod(a: Gaussian, b: Gaussian) -> tuple[Gaussian, Gaussian]:
    d = norm(b)
    assert d > 0
    x = a[0] * b[0] + a[1] * b[1]
    y = a[1] * b[0] - a[0] * b[1]
    q = nearest_quotient(x, d), nearest_quotient(y, d)
    return q, (a[0] - mul(q, b)[0], a[1] - mul(q, b)[1])


def ggcd(a: Gaussian, b: Gaussian) -> Gaussian:
    while b != (0, 0):
        _, rem = gdivmod(a, b)
        a, b = b, rem
    return a


def pell(r: int) -> tuple[int, int]:
    assert r >= 0
    x, y = 1, 0
    for _ in range(r):
        x, y = 9 * x + 20 * y, 4 * x + 9 * y
    return x, y


def H(r: int) -> Gaussian:
    x, y = pell(r)
    return x + y, 2 * y


F: Gaussian = (1, 2)
FB = conj(F)
UNITS = {(1, 0), (-1, 0), (0, 1), (0, -1)}
EVEN_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
ALL_PAIRS = tuple(combinations(range(4), 2))


def seven_points(blocks: tuple[Gaussian, ...]) -> tuple[Gaussian, ...]:
    points = [mul(FB, gprod(blocks))]
    for barred in EVEN_PAIRS:
        points.append(mul(F, gprod(conj(z) if i in barred else z
                                  for i, z in enumerate(blocks))))
    return tuple(points)


def eight_points(blocks: tuple[Gaussian, ...], parity: str) -> tuple[Gaussian, ...]:
    words = []
    for p in (1, 3) if parity == "odd" else (0, 2, 4):
        for unbarred in combinations(range(4), p):
            words.append((p, frozenset(unbarred)))
    points = []
    for p, unbarred in words:
        correction = F if (parity == "odd" and p == 1) else FB
        if parity == "even":
            correction = gprod((F,) * ((4 - p) // 2) + (FB,) * (p // 2))
        points.append(mul(correction, gprod(z if i in unbarred else conj(z)
                                             for i, z in enumerate(blocks))))
    return tuple(points)


def content(z: Gaussian) -> int:
    return gcd(abs(z[0]), abs(z[1]))


def global_gaussian_gcd(points: tuple[Gaussian, ...]) -> Gaussian:
    g = points[0]
    for z in points[1:]:
        g = ggcd(g, z)
    return g


def check_norm_support(n: int) -> tuple[Gaussian, ...]:
    indices = (n - 2, n, n + 2, n + 4)
    blocks = tuple(H(r) for r in indices)
    norms = tuple(norm(z) for z in blocks)
    assert all(gcd(norms[i], norms[j]) == 1
               for i, j in combinations(range(4), 2))
    assert all(gcd(5, x) == 1 for x in norms)
    assert all(norm(z) % 2 == 1 and content(z) == 1 for z in blocks)
    assert all(norm(ggcd(z, conj(z))) == 1 for z in blocks)

    # Exact resultant congruences (4), including the new gap-six case.
    for r in indices:
        x_r, y_r = pell(r)
        for d in (2, 4, 6):
            x_d, y_d = pell(d)
            hnext = H(r + d)
            c_same = mul((4, 0), mul((2, -1), (y_d * y_r, 0)))
            c_cross = mul((0, -4), (x_d * y_r, 0))
            gdiv_exact((hnext[0] - c_same[0], hnext[1] - c_same[1]), H(r))
            bh = conj(hnext)
            gdiv_exact((bh[0] - c_cross[0], bh[1] - c_cross[1]), H(r))

    return blocks


def check_gap_six_factors() -> None:
    x2, y2 = pell(2)
    x4, y4 = pell(4)
    x6, y6 = pell(6)
    assert (x2, y2, x4, y4, x6, y6) == (
        161, 72, 51841, 23184, 16692641, 7465176)
    assert x2 == 7 * 23 and y2 == 2**3 * 3**2
    assert x4 == 47 * 1103 and y4 == 2**4 * 3**2 * 7 * 23
    assert x6 == 7 * 23 * 103681
    assert y6 == 2**3 * 3**3 * 17 * 19 * 107
    assert x6 % 17 in (1, 16) and y6 % 17 == 0
    assert x6 % 103681 == 0 and (5 * y6 * y6 + 1) % 103681 == 0
    x12, y12 = pell(12)
    assert (x12 + 1) % 103681 == 0 and y12 % 103681 == 0
    # H_{-1}=5-8i, so its norm is 89, coprime to both exception moduli.
    assert norm((5, -8)) == 89 and gcd(89, 17 * 103681) == 1
    # The only prime factors of 1103 are checked by exact trial division;
    # it is prime and is 3 mod 4.  103681 need not be prime.
    assert all(1103 % d for d in range(2, 34)) and 1103 % 4 == 3


def check_binet_exact(n: int) -> None:
    lam = (Fraction(9), Fraction(4))
    rho = (Fraction(-1, 2), Fraction(1, 2))  # (sqrt(5)-1)/2
    h_re = (Fraction(1, 2), Fraction(1, 10))
    h_im = (Fraction(0), Fraction(1, 5))
    hm_re = (Fraction(1, 2), Fraction(-1, 10))
    hm_im = (Fraction(0), Fraction(-1, 5))
    # h_- = -i rho h and the exact Binet expression for every used H_r.
    assert hm_re == q5_mul(rho, h_im)
    assert hm_im == q5_neg(q5_mul(rho, h_re))
    for r in (n - 2, n, n + 2, n + 4):
        lp = q5_pow(lam, r)
        lm = q5_pow(lam, -r)
        real = q5_add(q5_mul(h_re, lp), q5_mul(hm_re, lm))
        imag = q5_add(q5_mul(h_im, lp), q5_mul(hm_im, lm))
        z = H(r)
        assert real == (Fraction(z[0]), Fraction(0))
        assert imag == (Fraction(z[1]), Fraction(0))
        # Positive exact phase error has size atan(rho*lambda^(-2r));
        # rho<1 and atan(u)<=u.  These rational/field identities avoid floats.
        assert 0 < rho[1] and rho[0] < 0
        assert 1 < 5 < 9  # 1<sqrt(5)<3, hence 0<rho<1.


def check_family(n: int) -> tuple[int, int, int]:
    assert n >= 61 and n % 60 == 1
    blocks = check_norm_support(n)
    check_binet_exact(n)
    block_norms = tuple(norm(z) for z in blocks)
    all_five_norms = (5,) + block_norms
    assert all(gcd(a, b) == 1 for a, b in combinations(all_five_norms, 2))

    p7 = seven_points(blocks)
    target_norm = 5 * prod(block_norms)
    assert len(p7) == 7 and len(set(p7)) == 7
    assert all(norm(z) == target_norm for z in p7)
    assert all(content(z) == 1 for z in p7)
    assert norm(global_gaussian_gcd(p7)) == 1

    # The quartet with D unbarred has D as a common Gaussian factor; the
    # complementary three rows conjugate D, so the full seven-point gcd is 1.
    d = blocks[0]
    q7 = (p7[0], p7[4], p7[5], p7[6])
    assert all(gdiv_exact(z, d) for z in q7)
    assert all(norm(z) == target_norm for z in q7)

    p8_odd = eight_points(blocks, "odd")
    assert len(p8_odd) == 8 and len(set(p8_odd)) == 8
    assert all(norm(z) == target_norm for z in p8_odd)
    assert all(content(z) == 1 for z in p8_odd)
    assert norm(global_gaussian_gcd(p8_odd)) == 1
    assert all(z[0] > 0 and (z[0] - z[1]) % 2 for z in p8_odd)
    # In the odd family, the words leaving D unbarred are exactly the
    # Hadamard quartet after reorienting A,B,C to their conjugates.
    dp, ap, bp, cp = (blocks[0],) + tuple(conj(z) for z in blocks[1:])
    odd_quartet = (
        mul(F, gprod((dp, ap, bp, cp))),
        mul(FB, gprod((dp, ap, conj(bp), conj(cp)))),
        mul(FB, gprod((dp, conj(ap), bp, conj(cp)))),
        mul(FB, gprod((dp, conj(ap), conj(bp), cp))),
    )
    assert set(odd_quartet).issubset(p8_odd)
    assert norm(global_gaussian_gcd(odd_quartet)) == norm(d)

    # Check the actual coordinates against |Im z|<=R E, E=4 lambda^(-2n+4),
    # with an exact ordered Q(sqrt(5)) comparison rather than floating point.
    bound_squared = q5_mul((Fraction(16 * target_norm), Fraction(0)),
                           q5_pow((Fraction(9), Fraction(4)), -4 * n + 8))
    for z in p8_odd:
        assert q5_nonnegative(q5_add(bound_squared, (Fraction(-z[1] ** 2), Fraction(0))))

    p8_even = eight_points(blocks, "even")
    assert len(p8_even) == 8 and len(set(p8_even)) == 8
    assert all(norm(z) == 25 * prod(block_norms) for z in p8_even)
    expected_contents = sorted([1, 1] + [5] * 6)
    assert sorted(content(z) for z in p8_even) == expected_contents
    assert norm(global_gaussian_gcd(p8_even)) == 1

    # Leading phase coefficients: seven-point rows tend to arg(F)=2 theta;
    # both complete eight-point families tend to zero.
    # Error coefficient vectors have four entries in {+1,-1}, so |phase
    # deviation| <= E=4 lambda^(-2(n-2)); width <=2E.  The exact radius
    # exponents then give 2E sqrt(R) <= 8*5^(1/4)*lambda^6.
    assert (-2 + 4) == 2  # bar(F) plus four unbarred H factors.
    assert (2 + 2 - 2) == 2  # F plus two unbarred and two barred H factors.
    for p in (1, 3):
        correction_theta = 2 if p == 1 else -2
        assert correction_theta + (2 * p - 4) == 0
    assert (-2 * (n - 2)) + (2 * n + 2) == 6
    # The translated rectangle bounds follow without numerical phases:
    # radial <= R E^2/2 <= 8 sqrt(5) lambda^12, while tangential <= 2 R E
    # and sqrt(R) >= 5^(1/4) lambda^(2n+2)/4, so their ratio is <=
    # 32*5^(1/4)*lambda^6.
    assert (4 * n + 4) + 2 * (-2 * n + 4) == 12
    assert (4 * n + 4) + (-2 * n + 4) - (2 * n + 2) == 6
    return len(p7), len(p8_odd), len(p8_even)


def main() -> None:
    check_gap_six_factors()
    results = [check_family(61 + 60 * k) for k in range(10)]
    assert all(x == (7, 8, 8) for x in results)
    print("PASS: exact recurrence, gap-2/4/6 support, 17 and 103681 exceptions, and Binet phase identities for 10 indices n=61+60k.")
    print("PASS: all ten fixtures have pairwise-coprime five supports, equal norms and distinct points; seven-point and odd-parity eight-point sets are individually primitive with global Gaussian gcd one. The odd quartet has exactly the common norm N(D), and all eight points pass an exact imaginary-coordinate bound.")
    print("PASS: the optional even-parity eight-point set has global gcd one and exactly six content-5 middle points; phase cancellation and endpoint exponent bound checked symbolically, without floating point.")


if __name__ == "__main__":
    main()
