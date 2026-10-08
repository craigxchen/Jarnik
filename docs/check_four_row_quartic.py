"""Exact u=1 check in Q(sqrt(2),sqrt(3)/2,i); standard library only."""

from fractions import Fraction as Q
from itertools import combinations


class K:
    """Basis a^e b^f i^g, with a^2=2, b^2=3/4, i^2=-1."""

    def __init__(self, x=0):
        self.c = x.c if isinstance(x, K) else (Q(x),) + (Q(0),) * 7

    @classmethod
    def coeffs(cls, c):
        x = cls()
        x.c = tuple(map(Q, c))
        return x

    def __add__(self, other):
        other = K(other)
        return K.coeffs(x + y for x, y in zip(self.c, other.c))

    __radd__ = __add__

    def __neg__(self):
        return K.coeffs(-x for x in self.c)

    def __sub__(self, other):
        return self + -K(other)

    def __rsub__(self, other):
        return K(other) + -self

    def __mul__(self, other):
        other = K(other)
        result = [Q(0)] * 8
        for j, x in enumerate(self.c):
            for k, y in enumerate(other.c):
                value = x * y
                for bit, square in [(1, Q(2)), (2, Q(3, 4)), (4, Q(-1))]:
                    if j & k & bit:
                        value *= square
                result[j ^ k] += value
        return K.coeffs(result)

    __rmul__ = __mul__

    def conjugate(self):
        return K.coeffs(-x if j & 4 else x for j, x in enumerate(self.c))

    def real(self):
        return self == self.conjugate()

    def __eq__(self, other):
        return self.c == K(other).c


def basis(j):
    return K.coeffs(int(k == j) for k in range(8))


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def multiply(p, q):
    r = [K(0)] * (len(p) + len(q) - 1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            r[j + k] = r[j + k] + a * b
    return trim(r)


def product_of_roots(roots):
    p = [K(1)]
    for root in roots:
        p = multiply(p, [-root, K(1)])
    return p


def monic_exact_division(p, q):
    assert q[-1] == 1
    p = list(p)
    result = [K(0)] * (len(p) - len(q) + 1)
    while len(p) >= len(q):
        shift, c = len(p) - len(q), p[-1]
        result[shift] = c
        for j, a in enumerate(q):
            p[shift + j] = p[shift + j] - c * a
        p = trim(p)
    assert p == [K(0)]
    return trim(result)


def main():
    a, b, imaginary = basis(1), basis(2), basis(4)
    triple = [a + imaginary, -a + imaginary,
              b - Q(3, 2) * imaginary, -b - Q(3, 2) * imaginary]
    private = [Q(1, 2) * a + 2 * imaginary, -Q(1, 2) * a + 2 * imaginary,
               3 * b - Q(1, 2) * imaginary, -3 * b - Q(1, 2) * imaginary]
    roots = triple + private
    assert len({r.c for r in roots}) == 8
    assert all(r != s.conjugate() for r in roots for s in roots)
    rows = [product_of_roots([private[j]] + [triple[k] for k in range(4) if k != j])
            for j in range(4)]
    assert rows[0] == [9 + Q(9, 2) * a * imaginary, 3 * a, K(3), Q(1, 2) * a, K(1)]
    assert rows[2] == [9 + 12 * b * imaginary, -2 * b, K(-2), -2 * b, K(1)]
    for j in (0, 2):
        assert rows[j + 1] == [((-1) ** k) * c.conjugate() for k, c in enumerate(rows[j])]
    for row in rows:
        assert all(c.real() for c in row[1:])
        assert not row[0].real()
    for j, k in combinations(range(4), 2):
        common = product_of_roots([triple[l] for l in range(4) if l not in (j, k)])
        norm_common = multiply(common, [c.conjugate() for c in common])
        pair = monic_exact_division(multiply(rows[j], [c.conjugate() for c in rows[k]]), norm_common)
        assert len(pair) == 5
        assert all(c.real() for c in pair[1:])
        assert not pair[0].real()
    assert len(product_of_roots(roots)) == 9
    print("Passed exact u=1 check: eight distinct conjugate-disjoint roots; four quartics; six constant-imaginary primitive pairs; lcm degree eight.")


if __name__ == "__main__":
    main()
