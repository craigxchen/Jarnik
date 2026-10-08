#!/usr/bin/env python3
"""Exact checks for positive_sum_full_cube_audit.md.

Only Python's standard library is used.  Gaussian integers are represented
by integer pairs, so all identities and inequalities are checked without
floating point arithmetic.
"""

from fractions import Fraction
from itertools import product
from math import gcd


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def mul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def conj(x):
    return (x[0], -x[1])


def norm(x):
    return x[0] * x[0] + x[1] * x[1]


def scale(x, n):
    return (n * x[0], n * x[1])


def primitive_conjugate_pair(k):
    a, b = k
    # gcd(a,b)=1 removes common split-prime factors; opposite parity
    # removes the common ramified factor 1+i.
    return (
        a != 0
        and b != 0
        and gcd(abs(a), abs(b)) == 1
        and (a - b) % 2 != 0
    )


def cube_points(d, kappas):
    rows = []
    for eps in product((1, -1), repeat=3):
        z = d
        for e, k in zip(eps, kappas):
            z = mul(z, k if e == 1 else conj(k))
        rows.append(z)
    return rows


def main():
    d = (2, 1)
    kappas = ((5, 2), (7, 2), (11, 4))
    assert all(primitive_conjugate_pair(k) for k in kappas)

    rows = cube_points(d, kappas)
    assert len(set(rows)) == 8
    rho2 = norm(d)
    for k in kappas:
        rho2 *= norm(k)
    assert all(norm(z) == rho2 for z in rows)

    s = (0, 0)
    for z in rows:
        s = add(s, z)
    expected = d
    for a, _ in kappas:
        expected = scale(expected, 2 * a)
    assert s == expected

    deficit = 8 * 8 * rho2 - norm(s)
    pair_deficit = 0
    for i in range(8):
        for j in range(i):
            pair_deficit += norm(add(rows[i], scale(rows[j], -1)))
    assert deficit == pair_deficit

    # Check (2) by exact cross multiplication for every block.
    for a, b in kappas:
        q2 = a * a + b * b
        assert deficit * q2 >= 64 * rho2 * b * b

    # Check the derived constants at C=1/2.  The certificate's conclusion
    # is rho <= (7/64)^3, while every nonzero Gaussian point has rho >= 1.
    bound = Fraction(7, 64) ** 3
    assert bound < 1
    assert Fraction(rho2, 1) > bound * bound

    # Power sums also factor, illustrating the exact scope of independence.
    for ell in (1, 2, 3, 4):
        lhs = (0, 0)
        for z in rows:
            p = (1, 0)
            for _ in range(ell):
                p = mul(p, z)
            lhs = add(lhs, p)
        rhs = (1, 0)
        dpow = (1, 0)
        for _ in range(ell):
            dpow = mul(dpow, d)
        rhs = dpow
        for k in kappas:
            kp = (1, 0)
            bp = (1, 0)
            for _ in range(ell):
                kp = mul(kp, k)
                bp = mul(bp, conj(k))
            rhs = mul(rhs, add(kp, bp))
        assert lhs == rhs

    print("positive-sum full-cube exact checks: PASS")
    print("rho^2 =", rho2, "sum deficit =", deficit, "C=1/2 bound < 1")


if __name__ == "__main__":
    main()
