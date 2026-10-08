"""Exact arithmetic fixtures for full_boolean_small_imaginary_circuits.md."""

from math import gcd
from fractions import Fraction
from itertools import product


def gaussian_power(a: int, b: int, exponent: int) -> tuple[int, int]:
    x, y = 1, 0
    for _ in range(exponent):
        x, y = x * a - y * b, x * b + y * a
    return x, y


def trim(p: list[Fraction]) -> list[Fraction]:
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    return trim([
        (p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
        for i in range(max(len(p), len(q)))
    ])


def sub(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    return add(p, [-x for x in q])


def mul(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)


def divide(p: list[Fraction], q: list[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    p = p[:]
    quotient = [Fraction(0)] * max(1, len(p) - len(q) + 1)
    while len(p) >= len(q) and p != [0]:
        shift = len(p) - len(q)
        coefficient = p[-1] / q[-1]
        quotient[shift] = coefficient
        p = sub(p, [Fraction(0)] * shift + [coefficient * x for x in q])
    return trim(quotient), p


def remainder(p: list[Fraction], q: list[Fraction]) -> list[Fraction]:
    return divide(p, q)[1]


def inverse(p: list[Fraction], modulus: list[Fraction]) -> list[Fraction]:
    r0, r1 = modulus, p
    s0, s1 = [Fraction(0)], [Fraction(1)]
    while r1 != [0]:
        quotient, _ = divide(r0, r1)
        r0, r1 = r1, sub(r0, mul(quotient, r1))
        s0, s1 = s1, sub(s0, mul(quotient, s1))
    assert len(r0) == 1
    return remainder([x / r0[0] for x in s0], modulus)


def check_pair_crt_leading_coefficients() -> None:
    """Canonical whole-block orientations force a degree-five shear."""
    F = Fraction
    u = [F(25), F(-6), F(1)]
    v = [F(25), F(6), F(1)]
    w = [F(25), F(0), F(1)]
    moduli = [u, v, w]
    modulus = mul(mul(u, v), w)
    roots = [[F(-3, 4), F(1, 4)],
             [F(3, 4), F(1, 4)],
             [F(0), F(1, 5)]]
    expected = [F(781, 45000), F(89, 5000),
                F(311, 18000), F(319, 18000),
                F(311, 18000), F(319, 18000),
                F(43, 2500), F(397, 22500)]

    def polynomial_crt(signs: tuple[int, int, int],
                       g: list[Fraction]) -> list[Fraction]:
        bases = [mul([F(2)], mul(g, w)),
                 mul([F(3)], mul(g, w)), [F(0)]]
        crt = [F(0)]
        for i, local_modulus in enumerate(moduli):
            rest, zero = divide(modulus, local_modulus)
            assert zero == [0]
            local_inverse = inverse(remainder(rest, local_modulus), local_modulus)
            residue = sub([signs[i] * x for x in roots[i]], bases[i])
            crt = add(crt, mul(residue, mul(rest, local_inverse)))
        crt = remainder(crt, modulus)
        for i, local_modulus in enumerate(moduli):
            residue = sub([signs[i] * x for x in roots[i]], bases[i])
            assert remainder(sub(crt, residue), local_modulus) == [0]
        return crt

    for signs, leading in zip(product((-1, 1), repeat=3), expected):
        crt = polynomial_crt(signs, [F(1), F(-1), F(1, 4)])
        assert len(crt) == 6 and crt[-1] == leading
        assert all(90000 % coefficient.denominator == 0 for coefficient in crt)
        c0 = polynomial_crt(signs, [F(0)])[-1]
        g0 = polynomial_crt(signs, [F(1)])[-1] - c0
        g1 = polynomial_crt(signs, [F(0), F(1)])[-1] - c0
        g2 = polynomial_crt(signs, [F(0), F(0), F(1)])[-1] - c0
        assert (g0, g1, g2) == (F(-1, 300), F(0), F(1, 12))
        assert leading == c0 - F(1, 300) + F(1, 48)
        assert 300 * c0 != int(300 * c0)


def check(h: int) -> None:
    a = 2 * (5**h + 1)
    g = 5 ** (2 * h)
    u = (a - 3) ** 2 + 16
    v = (a + 3) ** 2 + 16
    w = a * a + 25

    assert a % 2 == 0 and a % 5 != 0
    assert u + v == 2 * w
    assert all(n % 4 == 1 for n in (u, v, w, g))
    assert gcd(u, v) == gcd(u, w) == gcd(v, w) == 1
    assert gcd(g, u * v * w) == 1
    assert gcd(a - 3, 4) == gcd(a + 3, 4) == gcd(a, 5) == 1

    x = (3 * g * w, g * (6 * w - v), 2 * g * w)
    y = (1, 2, 1)
    d = lambda i, j: x[i] * y[j] - x[j] * y[i]
    assert min(x) > 0
    assert (d(0, 1), d(0, 2), d(1, 2)) == (g * v, g * w, g * u)
    assert y[0] * u - y[1] * w + y[2] * v == 0

    ga, gb = gaussian_power(2, 1, 2 * h)
    assert ga * ga + gb * gb == g and gcd(gb, g) == 1
    r = (ga * pow(gb, -1, g)) % g
    if r > g // 2:
        r -= g
    sheared_x = tuple(xi + r * yi for xi, yi in zip(x, y))
    assert min(sheared_x) > 0
    assert tuple(
        sheared_x[i] * y[j] - sheared_x[j] * y[i]
        for i, j in ((0, 1), (0, 2), (1, 2))
    ) == (g * v, g * w, g * u)
    for xi, yi in zip(sheared_x, y):
        # The real and imaginary coordinates after multiplying by
        # the conjugate of (2+i)^(2h) are both multiples of its norm.
        assert (xi * ga + yi * gb) % g == 0
        assert (yi * ga - xi * gb) % g == 0
        assert (xi * xi + yi * yi) % g == 0
        assert yi % 5 != 0

    if h == 1:
        assert a == 12 and v == 241
        assert (x[0] ** 2 + y[0] ** 2) % v == 206
        assert r == 7
        assert (sheared_x[0] ** 2 + y[0] ** 2) % v == 88


def main() -> None:
    for h in range(1, 13):
        check(h)
    check_pair_crt_leading_coefficients()
    print("PASS: 12 split-norm fixtures and 8 exact polynomial CRT coefficients")


if __name__ == "__main__":
    main()
