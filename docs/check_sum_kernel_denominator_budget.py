"""Exact sum-kernel denominator checks on actual Gaussian circles.

These finite tuples are not asserted to be endpoint families. The note
proves the denominator identity for arbitrary prime powers.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import atan2, factorial, gcd, lcm, prod

from check_reciprocal_chord_content import add, conj, divide, gp, gprod, mul, norm


def matchings(vertices):
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for second in vertices[1:]:
        for tail in matchings(tuple(i for i in vertices if i not in (first, second))):
            yield ((first, second),) + tail


MATCHINGS = list(matchings(tuple(range(8))))
assert len(MATCHINGS) == 105
PERMS = list(permutations(range(8)))
SIGNS = [(-1) ** sum(s[i] > s[j] for i in range(8) for j in range(i + 1, 8))
         for s in PERMS]
PRIMES = [(10, 1), (14, 1), (20, 1)]
NORMS = [norm(p) for p in PRIMES]

for exponents in ((1, 1, 1), (2, 1, 1), (1, 2, 2)):
    rows = [tuple(e * b for e, b in zip(exponents, bits))
            for bits in product(range(2), repeat=3)]
    points = [gprod(mul(gp(p, a), gp(conj(p), e - a))
                    for p, a, e in zip(PRIMES, row, exponents))
              for row in rows]
    angles = [atan2(z[1], z[0]) for z in points]
    assert max(angles) - min(angles) < 1
    assert len(set(points)) == 8 and len({norm(z) for z in points}) == 1
    x, t, q = {}, {}, {}
    for i in range(8):
        x[i, i], t[i, i], q[i, i] = 1, 0, 1
    for i, j in combinations(range(8), 2):
        core = gprod(mul(gp(p, min(rows[i][k], rows[j][k])),
                         gp(conj(p), e - max(rows[i][k], rows[j][k])))
                     for k, (p, e) in enumerate(zip(PRIMES, exponents)))
        real = divide(add(points[i], points[j]), mul((2, 0), core))
        imag = divide(add(points[i], tuple(-a for a in points[j])), mul((0, 2), core))
        assert real[1] == imag[1] == 0
        xx, tt = abs(real[0]), abs(imag[0])
        qq = prod(pn ** abs(rows[i][k] - rows[j][k]) for k, pn in enumerate(NORMS))
        assert xx > 0 and tt > 0 and xx * xx + tt * tt == qq
        assert gcd(xx, tt) == 1
        x[i, j] = x[j, i] = xx
        t[i, j] = t[j, i] = tt
        q[i, j] = q[j, i] = qq

    X = prod(x[i, j] for i, j in combinations(range(8), 2))
    T = prod(t[i, j] for i, j in combinations(range(8), 2))
    LX = lcm(*(prod(x[i, j] for i, j in matching) for matching in MATCHINGS))
    assert X % LX == 0 and T % (X // LX) == 0
    expected_lcd = LX * LX
    actual_lcd = 1
    permanent_numerator = 0
    squared_determinant_numerator = 0
    for sigma, sign in zip(PERMS, SIGNS):
        AA = prod(pn ** (sum(abs(rows[i][k] - rows[sigma[i]][k]) for i in range(8)) // 2)
                  for k, pn in enumerate(NORMS))
        XX = prod(x[i, sigma[i]] for i in range(8))
        assert AA * AA == prod(q[i, sigma[i]] for i in range(8))
        actual_lcd = lcm(actual_lcd, XX // gcd(AA, XX))
        assert expected_lcd * AA % XX == 0
        integer_term = expected_lcd * AA // XX
        assert integer_term >= expected_lcd
        permanent_numerator += integer_term
        squared_determinant_numerator += sign * integer_term * integer_term
    assert actual_lcd == expected_lcd
    assert actual_lcd % 2 == 0  # These cases exercise the shifted modulus at two.
    assert squared_determinant_numerator * X * X == T * T * expected_lcd * permanent_numerator

    permanent = Fraction(permanent_numerator, expected_lcd)
    assert permanent > factorial(8)
    taylor_terms = [Fraction(factorial(7) * t[i, j] ** 2, q[i, j])
                    for i, j in combinations(range(8), 2)]
    polynomial1 = factorial(8) + sum(taylor_terms)
    remainder1 = permanent - polynomial1
    H1 = lcm(expected_lcd, *(term.denominator for term in taylor_terms))
    assert remainder1 > 0 and (H1 * remainder1).denominator == 1

print("Three actual Gaussian eight-tuples: all40,320 permutation products checked each.")
print("Exact common denominator L_X^2 verified, including nontrivial two-adic factors.")
print("Normalized symmetric Cauchy/Borchardt identity and positive Taylor remainder pass.")
print("The test tuples are not asserted to be endpoint families.")
