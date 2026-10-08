"""Exact arithmetic audit of actual_endpoint_spectral_nonreal_difference.md."""

from fractions import Fraction as F
from math import gcd, isqrt


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for j, v in enumerate(a):
        for k, w in enumerate(b):
            out[j + k] += v * w
    return out


def rootpoly(roots):
    result = [F(1)]
    for root in roots:
        result = pmul(result, [-root, F(1)])
    return result


def reflect(p):
    return [v * (-1) ** j for j, v in enumerate(p)]


def at_i(p):
    result = (F(0), F(0))
    for v in reversed(p):
        result = gmul(result, (0, 1))
        result = (result[0] + v, result[1])
    return result


def check(count=32):
    x, y = 0, 1
    assert 1351 ** 2 - 75 * 156 ** 2 == 1
    for k in range(1, count + 1):
        x, y = 1351 * x + 156 * y, 11700 * x + 1351 * y
        assert y * y - 75 * x * x == 1
        assert x % 2 == 0 and y % 2 == 1 and gcd(x, y) == 1
        assert y % 5 == 1 and y % 13 == (-1) ** k % 13
        assert x >= 156 and 8 * x < y < 9 * x
        blocks = [(x, y - 8 * x), (5 * x, y), (x, y + 8 * x)]
        norms = [norm(z) for z in blocks]
        for z, n in zip(blocks, norms):
            assert z[0] > 0 and z[1] > 0 and gcd(*z) == 1 and n % 2
        for i in range(3):
            for j in range(i):
                assert gcd(norms[i], norms[j]) == 1
        assert norms == [140*x*x - 16*x*y + 1, 100*x*x + 1,
                         140*x*x + 16*x*y + 1]
        assert 10*x == isqrt(norms[1]) and norms[1] != (10*x)**2
        X = x * (200*x*x + 7)
        G = gmul(gmul(blocks[0], blocks[1]), blocks[2])
        assert G == (-X, -y)
        assert gcd(X, y) == 1 and norm(G) % 2
        assert norm(G) == norms[0] * norms[1] * norms[2]
        assert 75 * (200*x*x + 7) == 200*y*y + 325
        # Exact premises of C <= (27/20) / sqrt(x).
        assert y <= 9*x and X >= 200*x**3 and X+y <= 225*x**3
        assert 2 * F(9, 200) * 15 == F(27, 20)
        q = F(y, x)
        slopes = [q-8, q/5, q+8]
        assert 0 < slopes[0] < slopes[1] < slopes[2]
        pp = rootpoly(slopes)
        pm = rootpoly([-a for a in slopes])
        assert pm == [-v for v in reflect(pp)]
        assert pmul(pp, reflect(pp)) == pmul(pm, reflect(pm))
        shared = [F(-1)]
        for a in slopes:
            shared = pmul(shared, [-a*a, F(0), F(1)])
        assert pmul(pp, reflect(pp)) == shared
        e1 = sum(slopes)
        e3 = slopes[0] * slopes[1] * slopes[2]
        assert e1 == 11*q/5 and e3 == q*(11+F(1, x*x))/5
        difference = [a-b for a, b in zip(pp, pm)]
        assert difference == [-2*e3, 0, -2*e1, 0]
        squared_imaginary_root = e3/e1
        assert squared_imaginary_root == 1+F(1, 11*x*x) > 1
        assert -2*e3 + (-2*e1)*(-squared_imaginary_root) == 0
        assert at_i(pp) == (-F(y, 5*x**3), F(X, 5*x**3))
        assert at_i(pm) == (F(y, 5*x**3), F(X, 5*x**3))
        assert at_i(difference) == (-F(2*y, 5*x**3), F(0))
        # Logarithmic derivatives at the actual evaluation point.
        A = sum(a/(1+a*a) for a in slopes)
        B = sum(1/(1+a*a) for a in slopes)
        D = q**4-126*q*q+4225
        assert A == 2*q*(q*q-63)/D+5*q/(q*q+25)
        assert B == 2*(q*q+65)/D+25/(q*q+25)
        pp_derivative = [j*pp[j] for j in range(1, len(pp))]
        pp_at_i = at_i(pp)
        minus_derivative = tuple(-v for v in at_i(pp_derivative))
        assert gmul((A, B), pp_at_i) == minus_derivative
        # Every sign row on the three factors, testing the kernel identity.
        signs = [[1 if mask & (1 << c) else -1 for c in range(3)]
                 for mask in range(8)]
        values = [[(-s*a/(1+a*a), -1/(1+a*a))
                   for s, a in zip(row, slopes)] for row in signs]
        weights = [a*a/(1+a*a)**2 for a in slopes]
        for c in range(3):
            assert weights[c] == F(blocks[c][0]**2*blocks[c][1]**2,
                                   norms[c]**2)
            assert F(norms[c]-1, norms[c]**2) <= weights[c] <= F(1, 4)
        for i in range(8):
            for j in range(8):
                actual = (F(0), F(0))
                for c in range(3):
                    left = tuple(values[i][c][h]-values[0][c][h]
                                 for h in range(2))
                    right = (values[j][c][0]-values[0][c][0],
                             -values[j][c][1]+values[0][c][1])
                    term = gmul(left, right)
                    actual = tuple(actual[h]+term[h] for h in range(2))
                expected = sum((signs[i][c]-signs[0][c]) *
                               (signs[j][c]-signs[0][c]) * weights[c]
                               for c in range(3))
                assert actual == (expected, 0)
    assert F(2*(75-63), 75**2-126*75+4225)+F(5, 75+25) == F(11, 100)
    assert F(2*(75+65), 75**2-126*75+4225)+F(25, 75+25) == F(19, 20)
    print("PASS: %d exact Pell fixtures; primitive products and pairwise "
          "coprime rational norms; spectral products, nonreal difference "
          "roots, complex evaluations, logarithmic derivatives, and "
          "8-row centered kernel identities." % count)


if __name__ == "__main__":
    check()
