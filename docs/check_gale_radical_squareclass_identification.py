"""Exact checks of metric, full half-angle, and radical field squareclasses."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, lcm, prod


def rational_square(value):
    return value >= 0 and isqrt(value.numerator)**2 == value.numerator and \
        isqrt(value.denominator)**2 == value.denominator


def multiply(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conjugate(z):
    return z[0], -z[1]


def norm(z):
    return z[0]**2+z[1]**2


def divide(z, w):
    numerator = multiply(z, conjugate(w))
    return tuple(value/norm(w) for value in numerator)


def half_angle(phase):
    h = (F(0), F(1)) if phase == (F(-1), F(0)) else (1+phase[0], phase[1])
    denominator = lcm(*(value.denominator for value in h))
    integral = [int(value*denominator) for value in h]
    content = gcd(*map(abs, integral))
    h = tuple(F(value//content) for value in integral)
    assert divide(h, conjugate(h)) == phase
    return h


def rational_gcd(values):
    denominator = lcm(*(value.denominator for value in values))
    return F(reduce(gcd, (int(value*denominator) for value in values)), denominator)


def check_metric(form, nodes):
    a, b, c, r = map(F, form)
    assert a*c-b*b == r*r > 0
    vectors = [(F(1), F(t)) for t in nodes]
    values = [a*x*x+2*b*x*y+c*y*y for x, y in vectors]
    determinants = [prod(v[0]*w[1]-v[1]*w[0] for j, w in enumerate(vectors) if j != i)
                    for i, v in enumerate(vectors)]
    lengths = [q**5/d**2 for q, d in zip(values, determinants)]
    tau = rational_gcd(lengths)
    cs = [length/tau for length in lengths]
    assert all(value.denominator == 1 for value in cs)
    images = [(a*x+b*y, r*y) for x, y in vectors]
    for i, j in combinations(range(7), 2):
        beta = multiply(images[i], conjugate(images[j]))
        assert norm(beta) == a*a*values[i]*values[j]
        phase = divide(beta, conjugate(beta))
        h = half_angle(phase)
        assert rational_square(norm(h)/(values[i]*values[j]))
        assert rational_square(cs[i]*cs[j]/(values[i]*values[j]))
        assert divide(beta, h)[1] == 0


def check_degree_64():
    # Six independent relative norm classes, including a rational anchor.
    nodes = [0, 1, 2, 4, 6, 10, 14]
    norms = [1+t*t for t in nodes]
    assert norms == [1, 2, 5, 17, 37, 101, 197]
    assert all(all(p % d for d in range(2, isqrt(p)+1)) for p in norms[1:])
    derivatives = [prod(t-u for j, u in enumerate(nodes) if i != j)
                   for i, t in enumerate(nodes)]
    # Divide the radical by sqrt(c_0); its coefficients are rational times
    # 1,sqrt(2),...,sqrt(197). This square has all six singleton characters.
    coefficients = {0: F(1 if derivatives[0] > 0 else -1)}
    for i in range(1, 7):
        coefficients[1 << (i-1)] = F(norms[i]**2*abs(derivatives[0]), derivatives[i])
    square = {}
    for left, a in coefficients.items():
        for right, b in coefficients.items():
            intersection = left & right
            scalar = prod(norms[j+1] for j in range(6) if intersection >> j & 1)
            mask = left ^ right
            square[mask] = square.get(mask, F(0)) + a*b*scalar
    support = [mask for mask, value in square.items() if value]
    stabilizers = [g for g in range(64)
                   if all(bin(g & mask).count('1') % 2 == 0 for mask in support)]
    assert stabilizers == [0]


if __name__ == '__main__':
    forms = [(1, 0, 1, 1), (1, 2, 5, 1), (2, 1, 5, 3),
             (F(2, 3), F(1, 3), F(5, 3), 1)]
    for form in forms:
        check_metric(form, range(-3, 4))
        check_metric(form, [0, 1, 2, 4, 6, 10, 14])
    for phase, expected in [((F(1), F(0)), 1), ((F(-1), F(0)), 1),
                            ((F(0), F(1)), 2), ((F(0), F(-1)), 2)]:
        assert norm(half_angle(phase)) == expected
    check_degree_64()
    print('PASS: 168 exact metric/half-angle/Gale pair comparisons, all units, and degree-64 character example.')
