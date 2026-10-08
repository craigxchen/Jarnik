"""Exact arithmetic certificates; examples are not endpoint families."""
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, isqrt

from check_eight_point_pfaffian_content import EDGES, MATCHINGS, pair_data, prod


def lcm(a, b):
    return a // gcd(a, b) * b


def haf_cauchy(xs):
    if not xs:
        return Fraction(1)
    a = xs[0]
    return sum((haf_cauchy(xs[1:j] + xs[j + 1:]) / (1 + a * b)
                for j, b in enumerate(xs[1:], 1)), Fraction())


def valuation(n, p):
    n = abs(n)
    assert n
    answer = 0
    while n % p == 0:
        answer += 1
        n //= p
    return answer


def audit_tuple(xs):
    data = pair_data([(x, 1) for x in xs])
    assert all(x > 0 for x, _ in data.values())
    norms = {edge: x * x + t * t for edge, (x, t) in data.items()}
    norm_products = [prod(norms[e] for e in matching) for matching, _ in MATCHINGS]
    O = reduce(gcd, norm_products)
    x_products = [prod(data[e][0] for e in matching) for matching, _ in MATCHINGS]
    LX = reduce(lcm, x_products, 1)
    X = prod(data[e][0] for e in EDGES)
    T_signed = prod(data[e][1] for e in EDGES)
    BX = X // LX
    assert X % LX == 0 and T_signed % BX == 0
    terms = []
    pf_terms = []
    for ((matching, sign), qm, xm) in zip(MATCHINGS, norm_products, x_products):
        am = isqrt(qm // O)
        assert am * am * O == qm
        terms.append(Fraction(am, xm))
        pf_terms.append(sign * Fraction(am * prod(data[e][1] for e in matching), xm * xm))
    assert reduce(lcm, (term.denominator for term in terms), 1) == LX
    K = sum(term * LX for term in terms)
    assert K.denominator == 1 and K > 0
    K = K.numerator
    pf_integer = sum(pf_terms) * LX * LX
    assert pf_integer.denominator == 1
    assert pf_integer == (T_signed // BX) * K
    square = Fraction(O * K * K, LX * LX)
    direct = prod(x * x + 1 for x in xs) * haf_cauchy(tuple(xs)) ** 2
    assert square == direct and square > 105 ** 2
    return O, LX, K, abs(T_signed), square


P = (2297907909463, 113054407839984, 1602230240264968,
     7550950497889920, 10293635455067632)
SCALE = 721355647620465
FIXED = (38, 2, 22, 4, 10)


def evaluate(poly, x):
    return sum(a * x ** j for j, a in enumerate(poly))


def multiply_poly(a, b):
    result = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def main():
    numerator = [Fraction(0)] * 5
    for i in range(5):
        term = [haf_cauchy(FIXED[:i] + FIXED[i + 1:])]
        for j, a in enumerate(FIXED):
            if i != j:
                term = multiply_poly(term, [1, a])
        numerator = [a + b for a, b in zip(numerator, term)]
    assert tuple(a * SCALE for a in numerator) == P
    derivative = tuple((j + 1) * P[j + 1] for j in range(4))
    assert evaluate(P, 12) % 19 == 0
    assert evaluate(derivative, 12) % 19 == 15
    assert SCALE % 19 and all((1 + a * 12) % 19 for a in FIXED)
    for start, step in ((2, 2), (8, 4), (26, 6), (40, 10)):
        audit_tuple(tuple(start + step * j for j in range(8)))
    for xs, count in (((20, 18, 2, 28, 22, 6, 38, 4), 3),
                      ((20, 18, 2, 28, 22, 6, 4, 14), 4)):
        O, LX, K, T, square = audit_tuple(xs)
        data = pair_data([(x, 1) for x in xs])
        assert sum(x % 19 == 0 for x, _ in data.values()) == count
        assert O % 19 and K % 19 and T % 19
        assert valuation(square.denominator, 19) == 2 * valuation(LX, 19)

    root = 12
    modulus = 19
    examples = []
    for exponent in range(1, 9):
        if exponent > 1:
            previous = modulus // 19
            digit = (-evaluate(P, root) // previous) * pow(evaluate(derivative, root), -1, 19) % 19
            root += digit * previous
        assert evaluate(P, root) % modulus == 0
        xe = root if root % 2 == 0 else root + modulus
        be = -pow(20, -1, modulus) % modulus
        if be % 2:
            be += modulus
        if valuation(20 * be + 1, 19) != exponent:
            be += 2 * modulus
        assert xe >= 2 and be >= 2 and xe % 2 == be % 2 == 0
        assert valuation(20 * be + 1, 19) == exponent
        xs = (20, be) + FIXED + (xe,)
        assert len(set(x % 19 for x in xs)) == 8
        assert [(i, j) for i, j in EDGES if (1 + xs[i] * xs[j]) % 19 == 0] == [(0, 1)]
        h6 = haf_cauchy(FIXED + (xe,))
        assert h6.denominator % 19 and valuation(h6.numerator, 19) >= exponent
        O, LX, K, T, square = audit_tuple(xs)
        assert O % 19 and T % 19
        assert valuation(LX, 19) == exponent
        assert valuation(K, 19) >= exponent
        assert square.denominator % 19 and (square - 105 ** 2).denominator % 19
        examples.append((exponent, xe, be))
        modulus *= 19
    print("PASS: exact quartic Hafnian numerator and simple root modulo19.")
    print("PASS:14 actual Gaussian tuples, exact Hafnian LCD and Pfaffian factorization.")
    print("PASS:final denominator has no cancellation in the three- and four-pole-pair cases.")
    print("PASS:8 prime-power final-denominator cancellations with19 not dividing any old residue.")
    print("Hensel examples (e,x_e,b_e):", examples)


if __name__ == "__main__":
    main()
