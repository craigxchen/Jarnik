"""Exact finite certificate; asymptotic heights are proved in the companion note."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt


def add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    u, v = q
    if x == u and y == -v:
        return None
    slope = (v - y) / (u - x) if p != q else 3 * x * x / (2 * y)
    xx = slope * slope - x - u
    return xx, slope * (x - xx) - y


def multiple(n, point):
    if n < 0:
        return multiple(-n, (point[0], -point[1]))
    answer = None
    while n:
        if n & 1:
            answer = add(answer, point)
        point = add(point, point)
        n //= 2
    return answer


def homogeneous(point):
    return (F(0), F(1), F(0)) if point is None else (*point, F(1))


def determinant(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def main():
    point = (F(3), F(5))
    assert point[1] ** 2 == point[0] ** 3 - 2
    counts = {}
    for prime in (5, 7):
        counts[prime] = 1 + sum((y * y - x * x * x + 2) % prime == 0
                               for x in range(prime) for y in range(prime))
    assert counts == {5: 6, 7: 7}
    # Good reduction plus prime-to-p torsion injection forces trivial torsion.
    assert 1728 % 5 and 1728 % 7

    frame_indices = [0, 1, 3, 7]
    frame = [multiple(n, point) for n in frame_indices]
    for triple in combinations(frame, 3):
        assert determinant(*(homogeneous(p) for p in triple)) != 0
    boundary_indices = frame_indices + [-a - b for a, b in combinations(frame_indices, 2)]
    assert boundary_indices == [0, 1, 3, 7, -1, -3, -7, -4, -8, -10]
    assert len(set(boundary_indices)) == 10
    boundary_points = [multiple(n, point) for n in boundary_indices]
    assert len(set(boundary_points)) == 10
    for a, b in combinations(frame_indices, 2):
        third = multiple(-a - b, point)
        assert determinant(homogeneous(multiple(a, point)),
                           homogeneous(multiple(b, point)), homogeneous(third)) == 0

    # Restriction of -K is 9[O]-([O]+[P]+[3P]+[7P]), of degree5 and sum -11P.
    assert 9 - len(frame_indices) == 5
    s = -sum(frame_indices)
    assert s == -11 and sum(boundary_indices) == 2 * s
    for n in range(12, 101):
        canonical_l = 4 * n * n + (n - s) ** 2
        for k in boundary_indices:
            difference = F((n - k) ** 2) - F(canonical_l, 5)
            predicted = -2 * (F(k) - F(s, 5)) * n + k * k - F(s * s, 5)
            assert difference == predicted
        assert sum((n - k) ** 2 for k in boundary_indices) - 2 * canonical_l == 56

    def denominator_root(n):
        x, _ = multiple(n, point)
        root = isqrt(x.denominator)
        assert root * root == x.denominator
        return root

    fixed_bad_support = 6
    for difference in {abs(a - b) for a, b in combinations(boundary_indices, 2)}:
        fixed_bad_support *= denominator_root(difference)

    def strip_bad(value):
        while True:
            common = gcd(value, fixed_bad_support)
            if common == 1:
                return value
            value //= common

    for n in range(12, 26):
        contacts = [strip_bad(denominator_root(n - k)) for k in boundary_indices]
        assert all(gcd(a, b) == 1 for a, b in combinations(contacts, 2))

    print("Good-prime point counts:", counts)
    print("Four rational blow-up centers are in general position.")
    print("Ten distinct boundary indices:", boundary_indices)
    print("All six secant third-intersection identities pass exactly.")
    print("Canonical divisor degree5, Abel sum -11P, and all height-polynomial identities pass.")
    print("All ten good-prime denominator contacts are pairwise coprime for N=12,...,25.")
    print("No Gaussian full-profile or endpoint family is asserted.")


if __name__ == "__main__":
    main()
