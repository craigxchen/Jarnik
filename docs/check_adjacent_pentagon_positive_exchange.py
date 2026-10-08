#!/usr/bin/env python3
"""Exact checks for adjacent_pentagon_positive_exchange_cancellation.md."""

from fractions import Fraction
from itertools import product
from math import gcd

from check_positive_pentagon_primitive_content import (
    gaussian_div, group_gcd, group_norm_gcd, integer_content,
    integer_valuation, norm, pentagon, product_ratio, sub, valuation,
)


def exchange(x, y, z):
    assert 0 < x < 1 and 0 < y < 1 and 0 < z < 1
    assert x * z + y > 1
    a, b = x.numerator, x.denominator
    c, d = y.numerator, y.denominator
    e, f = z.numerator, z.denominator
    EL, ER = a * d + b * c - b * d, c * f + d * e - d * f
    T = a * d * e + b * c * f - b * d * f
    B0 = d * (d - c) * (b - a) * (f - e)
    M = EL * ER
    assert EL > 0 and ER > 0 and T > 0
    assert M == B0 + c * T
    _, factorsL, FL = pentagon(a, b, c, d)
    _, factorsR, FR = pentagon(c, d, e, f)
    assert EL == factorsL[0] * factorsL[1] * factorsL[2] * FL
    assert ER == factorsR[0] * factorsR[1] * factorsR[2] * FR
    qL, qR = Fraction(EL, a * d), Fraction(ER, d * e)
    W, r = qL * qR, Fraction(B0, M)
    A, B = W * r, W * (1 - r)
    assert W == A + B and A > 0 and B > 0
    P, Q, k = W.numerator, W.denominator, r.denominator
    J = gcd(A.numerator, B.numerator)
    assert J == P // gcd(P, k)
    U, S = gcd(M, a * e * d * d), gcd(M, B0)
    assert P == M // U and k == M // S
    assert J == S // gcd(S, U)
    G = gcd(qL.numerator, qR.denominator) * gcd(qR.numerator, qL.denominator)
    assert P == qL.numerator * qR.numerator // G
    assert Q == qL.denominator * qR.denominator // G
    s = (1 - x) / (x + y - 1)
    rr = (z + y - 1) / ((1 - y) * (1 - z))
    assert rr > s > 0
    return FL, FR, qL, qR, W, r, A, B, J


def gap(t):
    return max(0, min(t[0], t[4]) - max(t[1:4]),
               min(t[1:4]) - max(t[0], t[4]))


def fixture(points, primes):
    N = norm(points[0])
    assert all(norm(w) == N for w in points)
    assert group_norm_gcd(points) == 1
    for i in range(6):
        for j in range(i + 1, 6):
            u, v = points[i], points[j]
            assert u[0] * v[1] - u[1] * v[0] > 0
            assert u[0] * v[0] + u[1] * v[1] > 0
    xyz = [product_ratio(points, [(i, i + 3), (i + 1, i + 2)],
                         [(i, i + 2), (i + 1, i + 3)]) for i in range(3)]
    FL, FR, qL, qR, W, r, A, B, J = exchange(*xyz)
    assert qL == product_ratio(points, [(0, 4), (2, 3)], [(0, 3), (2, 4)])
    assert qR == product_ratio(points, [(1, 5), (2, 3)], [(1, 3), (2, 5)])
    assert r == product_ratio(points, [(0, 1), (4, 5)], [(0, 4), (1, 5)])
    shared_gcd = group_gcd(points[2:4])
    shared_content = integer_content([gaussian_div(sub(points[3], points[2]), shared_gcd)])
    checks = []
    for prime in primes:
        p = norm(prime)
        t = [valuation(w, prime) for w in points]
        hL, hR = gap(t[:5]), gap(t[1:])
        if not hL or not hR:
            continue
        def rho(i, j):
            return valuation(sub(points[i], points[j]), prime) - t[i] if t[i] == t[j] else 0
        r04, r15, r23 = rho(0, 4), rho(1, 5), rho(2, 3)
        width = abs(t[2] - t[3])
        assert integer_valuation(qL.numerator, p) == hL + r04 + r23
        assert integer_valuation(qR.numerator, p) == hR + r15 + r23
        assert qL.denominator % p and qR.denominator % p
        assert integer_valuation(W.numerator, p) == hL + hR + r04 + r15 + 2 * r23
        assert integer_valuation(r.denominator, p) == hL + hR + width + r04 + r15
        assert integer_valuation(J, p) == 2 * r23
        assert integer_valuation(shared_content, p) == r23
        assert integer_valuation(FL * FR, p) == hL + hR + r04 + r15
        checks.append((p, t, hL, hR, width, r23))
    assert checks
    assert 4 * r.denominator ** 2 <= norm(sub(points[0], points[4])) * norm(sub(points[1], points[5]))
    assert 16 * norm(sub(points[5], points[0])) ** 2 > N
    print(f"fixture N={N}: F=({FL},{FR}), W={W}, r={r}, A={A}, B={B}, J={J}; {checks}")


def main():
    fractions = sorted({Fraction(a, b) for b in range(2, 15) for a in range(1, b)})
    checked = 0
    for x, y, z in product(fractions, repeat=3):
        if x * z + y > 1:
            exchange(x, y, z)
            checked += 1
    fixture(((-107, -54), (-98, -69), (-91, -78), (-78, -91), (-69, -98), (-54, -107)),
            ((2, 1), (3, 2), (4, 1)))
    fixture(((323, 36), (312, 91), (300, 125), (280, 165), (260, 195), (253, 204)),
            ((2, 1), (3, 2)))
    fixture(((1170, 65), (1167, 106), (1150, 225), (650, 975), (538, 1041), (510, 1055)),
            ((2, 1), (3, 2)))
    print(f"PASS: {checked} positive rational six-point systems and three actual ordered primitive fixtures")


if __name__ == "__main__":
    main()
