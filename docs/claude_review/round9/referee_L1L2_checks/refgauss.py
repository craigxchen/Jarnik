"""Referee's own exact Gaussian-integer toolkit (independent of the author's glib.py).

Gaussian integers are pairs (a, b) of Python ints meaning a + b i.
"""
from math import atan2, pi, isqrt, log


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def gsub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def gadd(z, w):
    return (z[0] + w[0], z[1] + w[1])


def _rdiv(a, b):
    # nearest integer to a/b for integers, b > 0
    return (2 * a + b) // (2 * b)


def gdivmod(z, w):
    n = gnorm(w)
    num = gmul(z, gconj(w))
    q = (_rdiv(num[0], n), _rdiv(num[1], n))
    r = gsub(z, gmul(q, w))
    return q, r


def gexactdiv(z, w):
    """z / w, asserting exactness."""
    n = gnorm(w)
    num = gmul(z, gconj(w))
    assert num[0] % n == 0 and num[1] % n == 0, (z, w)
    return (num[0] // n, num[1] // n)


def gdivides(w, z):
    n = gnorm(w)
    num = gmul(z, gconj(w))
    return num[0] % n == 0 and num[1] % n == 0


def ggcd(z, w):
    while w != (0, 0):
        _, r = gdivmod(z, w)
        z, w = w, r
    return z


def gpow(z, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, z)
    return r


UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # i^0..i^3


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def gauss_prime(p):
    """A Gaussian prime of norm p (p = 1 mod 4)."""
    assert p % 4 == 1
    c = 2
    while pow(c, (p - 1) // 2, p) != p - 1:
        c += 1
    x = pow(c, (p - 1) // 4, p)
    g = ggcd((p, 0), (x, 1))
    assert gnorm(g) == p
    return g


def all_points(N):
    """All lattice points on x^2 + y^2 = N, by brute force (for moderate N)."""
    pts = []
    r = isqrt(N)
    for x in range(-r, r + 1):
        y2 = N - x * x
        y = isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    return pts


def windows(N, M, Cmax):
    """All windows of M angularly consecutive lattice points with arc constant <= Cmax.
    Arc constant = (theta_last - theta_first) * N^(1/4)."""
    pts = all_points(N)
    pts.sort(key=lambda z: atan2(z[1], z[0]) % (2 * pi))
    n = len(pts)
    ang = [atan2(z[1], z[0]) % (2 * pi) for z in pts]
    s = N ** 0.25
    out = []
    for i in range(n):
        j = (i + M - 1) % n
        d = (ang[j] - ang[i]) % (2 * pi)
        C = d * s
        if C <= Cmax:
            out.append(([pts[(i + t) % n] for t in range(M)], C))
    return out


def normalise(cl):
    G = cl[0]
    for z in cl[1:]:
        G = ggcd(G, z)
    u = [gexactdiv(z, G) for z in cl]
    return u, G


def vpi(z, p_):
    """valuation of Gaussian integer z at Gaussian prime p_."""
    v = 0
    while z != (0, 0) and gdivides(p_, z):
        z = gexactdiv(z, p_)
        v += 1
    return v
