"""Exact checks for the five-point affine-relation lattice example."""

from fractions import Fraction as Q
from functools import reduce
from itertools import combinations
from math import gcd, isqrt

FACTORS = [
    (5, -3, 4), (3085, -2117, -2244), (65, -56, 33),
    (145, 144, -17), (364033085, 59073189, 122657188),
    (1885, -1637, -2244), (984115, 2144984, 43923),
    (1483885, -1069488, -861101),
]


def gadd(z, w):
    return z[0] + w[0], z[1] + w[1]


def gsub(z, w):
    return z[0] - w[0], z[1] - w[1]


def gmul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gconj(z):
    return z[0], -z[1]


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def gdiv_exact(z, w):
    n = gnorm(w)
    q = gmul(z, gconj(w))
    assert q[0] % n == 0 and q[1] % n == 0
    return q[0] // n, q[1] // n


def ground_div(z, w):
    """Nearest quotient in Z[i], with exact integer tie handling."""
    n = gnorm(w)
    q = gmul(z, gconj(w))

    def nearest(a):
        lo = a // n
        return lo + (2 * (a - lo * n) >= n)

    return nearest(q[0]), nearest(q[1])


def ggcd(z, w):
    while w != (0, 0):
        q = ground_div(z, w)
        z, w = w, gsub(z, gmul(q, w))
    return z


def primitive_points(n):
    qs = [(d * n - a, -b) for d, a, b in FACTORS]
    product = reduce(gmul, qs, (1, 0))
    points = [gconj(product)]
    for j in range(4):
        h = reduce(gmul, [qs[4 + j]] + [qs[k] for k in range(4) if k != j], (1, 0))
        points.append(gdiv_exact(gmul(gconj(product), h), gconj(h)))
    common = reduce(ggcd, points)
    return [gdiv_exact(z, common) for z in points]


def determinant(u, v):
    return u[0] * v[1] - u[1] * v[0]


def pmul(f, h):
    ans = [(0, 0)] * (len(f) + len(h) - 1)
    for i, z in enumerate(f):
        for j, w in enumerate(h):
            ans[i + j] = gadd(ans[i + j], gmul(z, w))
    return ans


def psub(f, h):
    size = max(len(f), len(h))
    return [gsub(f[i] if i < len(f) else (0, 0),
                 h[i] if i < len(h) else (0, 0)) for i in range(size)]


def coordinate_polynomial(points, j):
    return [z[j] for z in points]


def determinant_polynomial(u, v):
    real_u, imag_u = coordinate_polynomial(u, 0), coordinate_polynomial(u, 1)
    real_v, imag_v = coordinate_polynomial(v, 0), coordinate_polynomial(v, 1)

    def integer_mul(f, h):
        ans = [0] * (len(f) + len(h) - 1)
        for i, a in enumerate(f):
            for j, b in enumerate(h):
                ans[i + j] += a * b
        return ans

    a, b = integer_mul(real_u, imag_v), integer_mul(imag_u, real_v)
    return [x - y for x, y in zip(a, b)]


def point_polynomials():
    factors = [[(-a, -b), (d, 0)] for d, a, b in FACTORS]
    conjugates = [[gconj(z) for z in f] for f in factors]
    anchor = reduce(pmul, conjugates, [(1, 0)])
    points = [anchor]
    for j in range(4):
        support = {4 + j} | {k for k in range(4) if k != j}
        oriented = [factors[k] if k in support else conjugates[k] for k in range(8)]
        points.append(reduce(pmul, oriented, [(1, 0)]))
    return points


def minor_polynomials():
    points = point_polynomials()
    ans = {}
    for i, j, k in combinations(range(5), 3):
        ans[(i, j, k)] = determinant_polynomial(psub(points[j], points[i]),
                                                 psub(points[k], points[i]))
    return ans


def qtrim(f):
    f = [Q(x) for x in f]
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def qremainder(f, h):
    f, h = qtrim(f), qtrim(h)
    while len(f) >= len(h) and f != [0]:
        c, shift = f[-1] / h[-1], len(f) - len(h)
        for j, a in enumerate(h):
            f[j + shift] -= c * a
        f = qtrim(f)
    return f


def polynomial_gcd(f, h):
    f, h = qtrim(f), qtrim(h)
    while h != [0]:
        f, h = h, qremainder(f, h)
    return [a / f[-1] for a in f]


def minors(points):
    ans = {}
    for i, j, k in combinations(range(len(points)), 3):
        ans[(i, j, k)] = determinant(gsub(points[j], points[i]), gsub(points[k], points[i]))
    return ans


def main():
    polys = {key: qtrim(coeffs) for key, coeffs in minor_polynomials().items()}
    degrees = {key: len(poly) - 1 for key, poly in polys.items()}
    assert all(degree == 4 for degree in degrees.values())
    common_poly = reduce(polynomial_gcd, polys.values())
    assert len(common_poly) == 1
    print("minor degrees", degrees, "monic polynomial gcd", common_poly)
    for n in (2, 10, 50, 100, 250):
        points = primitive_points(n)
        radius2 = gnorm(points[0])
        assert all(gnorm(z) == radius2 for z in points)
        ds = minors(points)
        content = reduce(gcd, (abs(x) for x in ds.values()))
        covol2_num = sum(x * x for x in ds.values())
        print(n, "radius~", isqrt(radius2), "minor-content", content,
              "covol~", isqrt(covol2_num) // content,
              "max-minor", max(abs(x) for x in ds.values()))


if __name__ == "__main__":
    main()
