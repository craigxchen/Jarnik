"""Exact CRT obstruction to a common low-degree scalar discriminant."""
from fractions import Fraction
from check_polynomial_real_jet_differential_divisor_bound import (
    add, scale, mul, sq, derivative, degree, remainder,
)


def divmodp(a, b):
    a, b = tuple(map(Fraction, a)), tuple(map(Fraction, b))
    q = (0,)
    while a != (0,) and len(a) >= len(b):
        term = (0,) * (len(a) - len(b)) + (a[-1] / b[-1],)
        q = add(q, term)
        a = add(a, scale(mul(term, b), -1))
    return q, a


def inverse(a, b):
    r0, r1, s0, s1 = a, b, (1,), (0,)
    while r1 != (0,):
        q, r = divmodp(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, add(s0, scale(mul(q, s1), -1))
    assert degree(r0) == 0
    return scale(s0, 1 / r0[0])


def substitute_square(a):
    result = [0] * (2 * len(a) - 1)
    result[::2] = a
    return tuple(result)


def main():
    data = [(404, 720, 85, 1164), (428, -720, 61, 396),
            (464, -720, 25, 324), (484, -720, 5, -276),
            (524, 720, -35, -1236)]
    P = (1,)
    for a in range(1, 7):
        P = mul(P, (a*a, 0, 1))
    B, M = (Fraction(0),), (Fraction(1),)
    rows = []
    for c, e, d, f in data:
        x, y = (e, 0, -c, 0, -45, 0, -1), (0, -f, 0, -d, 0, -1)
        rows.append((x, y))
        assert add(sq(x), sq(y)) == P
        m = (f, d, 1)
        assert f != 0
        assert all(f-d*a*a+a**4 != 0 for a in range(1, 7))
        r = mul((0, 4), sq((c, 90, 3)))
        correction = remainder(mul(add(r, scale(B, -1)), inverse(M, m)), m)
        B, M = add(B, mul(M, correction)), mul(M, m)
    assert degree(B) == 9 and B[-1] == Fraction(389, 13440000)
    assert B[0] == Fraction(370384723467, 4375)
    assert M[0] == 50947247932416
    C = add(B, scale(M, -B[0] / M[0]))
    assert C[0] == 0 and degree(C) == 10
    assert C[-1] == Fraction(-67, 40320000)
    CT = substitute_square(C)
    K = add(scale(mul(CT, P), 4), scale(sq(derivative(P)), -1))
    for x, y in rows:
        assert remainder(K, y) == (0,)
        assert degree(add(x, scale(rows[0][0], -1))) <= 2
        assert degree(add(y, scale(rows[0][1], -1))) <= 3
    # Factor supports certify polynomial gcd one: each distinct Gaussian
    # linear factor is absent from at least one row.
    subsets = [{4, 6}, {1, 3, 6}, {1, 4, 5}, {2, 3, 5}, {1, 2, 3, 4}]
    for S, (x, y) in zip(subsets, rows):
        real, imag = (1,), (0,)
        for a in range(1, 7):
            factor_real, factor_imag = (a,), (0, 1 if a in S else -1)
            real, imag = (add(mul(real, factor_real), scale(mul(imag, factor_imag), -1)),
                          add(mul(real, factor_imag), mul(imag, factor_real)))
        sign = (-1) ** len(S)
        assert scale(real, sign) == x and scale(imag, sign) == y
    for a in range(1, 7):
        assert any(a in S for S in subsets) and any(a not in S for S in subsets)
    print('Verified: primitive five-row fixture; pairwise-coprime CRT moduli; '
          'unique common scalar remainder has degree 20 and leading coefficient '
          '-67/40320000.')


if __name__ == '__main__':
    main()
