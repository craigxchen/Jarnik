#!/usr/bin/env python3
"""Small exact audit of the dyadic Runge truncation."""
from fractions import Fraction
from itertools import combinations


def mul(a, b, lim):
    out = [Fraction(0) for _ in range(min(lim, len(a) + len(b) - 1))]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < lim:
                out[i + j] += x * y
    return out


def trunc_sqrt(cs, k):
    f = [Fraction(1)]
    for c in cs:
        g = [Fraction(0)] * (k + 1)
        for r in range(k + 1):
            # (1-t)^(1/2) coefficient, with t=cX^{-1}
            x = Fraction(1)
            for i in range(r):
                x *= Fraction(1, 2) - i
            import math
            x /= math.factorial(r)
            g[r] = x * (-c)**r
        f = mul(f, g, k + 1)
    # f[r] multiplies X^(k-r)
    return list(reversed(f))


def poly_mul(a, b):
    z = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): z[i + j] += x * y
    return z


def main():
    cases = 0
    for k in range(1, 5):
        kcases = 0
        for D in range(2 * k, 13):
            for cs in combinations(range(1, D + 1), 2 * k):
                if kcases >= 125: break
                s = trunc_sqrt(cs, k)
                q = 2 ** (2 * k - 1)
                assert all((q * x).denominator == 1 for x in s)
                for a, x in enumerate(reversed(s)):
                    assert abs(x) <= (2 * k * D) ** a
                p = [Fraction(1)]
                for c in cs: p = poly_mul(p, [Fraction(-c), Fraction(1)])
                ss = poly_mul(s, s)
                r = [q*q*x for x in p]
                r[:len(ss)] = [r[i] - q*q*ss[i] for i in range(len(ss))]
                assert all(x.denominator == 1 for x in r)
                assert all(abs(x) <= (k + 2) * (8 * k * D) ** (2 * k) for x in r)
                # leading terms through degree k cancel; generic distinct roots leave R nonzero
                assert all(r[i] == 0 for i in range(k, len(r)))
                assert any(r[i] for i in range(k))
                cases += 1
                kcases += 1
            if kcases >= 125: break
    # Repeated-root control: P is a square, so the cleared R really vanishes.
    s = trunc_sqrt((5, 4, 5, 4), 2)
    p = poly_mul(poly_mul([Fraction(-5), Fraction(1)], [Fraction(-4), Fraction(1)]),
                 poly_mul([Fraction(-5), Fraction(1)], [Fraction(-4), Fraction(1)]))
    assert poly_mul(s, s) == p
    print(f"PASS: {cases} distinct-root truncations; repeated-root square control")


if __name__ == '__main__': main()
