"""Exact rational four-row quartic certificate; standard library only."""

from fractions import Fraction as Q
from itertools import combinations


class G:
    def __init__(self, real=0, imag=0):
        self.a, self.b = (real.a, real.b) if isinstance(real, G) else (Q(real), Q(imag))

    def __add__(self, other):
        other = G(other)
        return G(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.a, -self.b)

    def __sub__(self, other):
        return self + -G(other)

    def __rsub__(self, other):
        return G(other) + -self

    def __mul__(self, other):
        other = G(other)
        return G(self.a * other.a - self.b * other.b,
                 self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def conjugate(self):
        return G(self.a, -self.b)

    def __eq__(self, other):
        other = G(other)
        return (self.a, self.b) == (other.a, other.b)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def multiply(p, q):
    result = [G(0)] * (len(p) + len(q) - 1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            result[j + k] = result[j + k] + a * b
    return trim(result)


def polynomial(roots):
    result = [G(1)]
    for root in roots:
        result = multiply(result, [-root, G(1)])
    return result


def monic_exact_division(p, q):
    assert q[-1] == 1
    p = list(p)
    result = [G(0)] * (len(p) - len(q) + 1)
    while len(p) >= len(q):
        shift, c = len(p) - len(q), p[-1]
        result[shift] = c
        for j, a in enumerate(q):
            p[shift + j] = p[shift + j] - c * a
        p = trim(p)
    assert p == [G(0)]
    return trim(result)


def main():
    parameters = [Q(2), -Q(51, 22), Q(11, 3), -Q(1, 17)]
    assert sum(parameters[j] * parameters[k] for j, k in combinations(range(4), 2)) == -6
    prod = Q(1)
    for r in parameters:
        prod *= r
    assert prod == 1
    triple = [G((1-r*r)/(1+r*r), 2*r/(1+r*r)) for r in parameters]
    B = -sum(z.b for z in triple)
    C = -sum(z.a * z.b for z in triple)
    assert B == -Q(107712, 232609)
    assert C == Q(28929719808, 54106946881)
    private = [G((z.a * z.b + C)/(z.b + B), z.b + B) for z in triple]
    roots = triple + private
    integer_factors = [
        (5, -3, 4), (3085, -2117, -2244), (65, -56, 33),
        (145, 144, -17), (364033085, 59073189, 122657188),
        (1885, -1637, -2244), (984115, 2144984, 43923),
        (1483885, -1069488, -861101),
    ]
    for z, (den, real, imag) in zip(roots, integer_factors):
        assert z == G(Q(real, den), Q(imag, den))
        assert den % 2 == 1 and (real - imag) % 2 == 1
    assert len({(z.a, z.b) for z in roots}) == 8
    assert all(z.b != 0 for z in roots)
    assert all(z != w.conjugate() for z in roots for w in roots)
    row_roots = [[private[j]] + [triple[k] for k in range(4) if k != j] for j in range(4)]
    rows = [polynomial(rr) for rr in row_roots]
    imaginary_constants = [
        Q(16334372208, 45504135625), Q(206448, 235625),
        -Q(5077968, 11183125), Q(5592048, 185485625),
    ]
    for row, tau in zip(rows, imaginary_constants):
        assert len(row) == 5 and row[-1] == 1
        assert [c.b for c in row] == [tau, Q(0), Q(0), Q(0), Q(0)]
    assert len(set(imaginary_constants)) == 4
    for j, k in combinations(range(4), 2):
        common_roots = [triple[l] for l in range(4) if l not in (j, k)]
        assert set((z.a, z.b) for z in row_roots[j]) & set((z.a, z.b) for z in row_roots[k]) == set((z.a, z.b) for z in common_roots)
        common = polynomial(common_roots)
        norm_common = multiply(common, [c.conjugate() for c in common])
        pair = monic_exact_division(multiply(rows[j], [c.conjugate() for c in rows[k]]), norm_common)
        assert len(pair) == 5
        assert [c.b for c in pair] == [imaginary_constants[j] - imaginary_constants[k], Q(0), Q(0), Q(0), Q(0)]
    L = polynomial(roots)
    assert len(L) == 9
    anchor = [c.conjugate() for c in L]
    for row, tau in zip(rows, imaginary_constants):
        quotient = monic_exact_division(anchor, [c.conjugate() for c in row])
        point = multiply(quotient, row)
        difference = trim([a-b for a, b in zip(point, anchor)])
        assert len(difference) == 5
        assert difference == [G(0, 2*tau) * c for c in quotient]
    print("Passed: four rational constant-imaginary quartics; six exact primitive pairs; eight independent roots; lcm degree eight; five circle polynomials with differences of degree four.")


if __name__ == "__main__":
    main()
