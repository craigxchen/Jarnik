"""Independent exact library for the integer-cotangent dictionary.

Written from scratch for the growth/cotangent.md write-up (does not import the
research-notes checkers).  Gaussian integers are pairs (a, b) = a + b i.
Everything is exact integer / Fraction arithmetic.
"""
from fractions import Fraction
from math import gcd, isqrt
from itertools import combinations
import random


# ---------------------------------------------------------------- Gaussian ops
def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def gsub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def gdivmod_round(z, w):
    """Nearest-integer Gaussian division z = q w + r with N(r) < N(w)."""
    n = gnorm(w)
    num = gmul(z, gconj(w))
    def rdiv(a, b):  # round a/b to nearest integer, b > 0
        return (2 * a + b) // (2 * b)
    q = (rdiv(num[0], n), rdiv(num[1], n))
    r = gsub(z, gmul(q, w))
    return q, r


def gexact(z, w):
    """Exact Gaussian division; asserts divisibility."""
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
        _, r = gdivmod_round(z, w)
        z, w = w, r
    return z


def ggcd_all(zs):
    g = (0, 0)
    for z in zs:
        g = ggcd(g, z)
    return g


def glcm(z, w):
    g = ggcd(z, w)
    return gexact(gmul(z, w), g)


def glcm_all(zs):
    out = (1, 0)
    for z in zs:
        out = glcm(out, z)
    return out


def gpow(z, e):
    out = (1, 0)
    for _ in range(e):
        out = gmul(out, z)
    return out


# ------------------------------------------------------------ rational helpers
def lcm(a, b):
    return a // gcd(a, b) * b if a and b else 0


def lcm_all(xs):
    out = 1
    for x in xs:
        out = lcm(out, x)
    return out


def factor(n):
    """Trial-division factorisation (small n only)."""
    n = abs(n)
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def vp(n, p):
    if n == 0:
        return 10 ** 9
    n = abs(n)
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def split_prime_pi(p):
    """A Gaussian prime above a split rational prime p = 1 mod 4."""
    assert p % 4 == 1
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b * b == b2:
            return (a, b)
    raise ValueError(p)


# ------------------------------------------------------- the cotangent chart
def edge_norm(Q, L):
    """n_e = (a^2+b^2)/eps with a=Q/g, b=L/g, g=gcd(Q,L)."""
    g = gcd(Q, L)
    a, b = Q // g, L // g
    eps = 2 if (a % 2 and b % 2) else 1
    assert (a * a + b * b) % eps == 0
    return (a * a + b * b) // eps


def pair_cot(x, y, L):
    num = x * y + L * L
    assert num % (x - y) == 0, (x, y, L)
    return num // (x - y)


def is_clique(xs, L):
    return all((x * y + L * L) % (x - y) == 0 for x, y in combinations(xs, 2))


def all_edge_lcm(xs, L):
    edges = list(xs) + [pair_cot(x, y, L) for x, y in combinations(xs, 2)]
    return lcm_all(edge_norm(q, L) for q in edges)


def primitive_tuple(xs, L):
    """Primitive Gaussian realisation of q_0=1, q_i=(x+iL)/(x-iL)."""
    dens = []
    nums = []
    for x in xs:
        h = (x, L)
        c = ggcd(h, gconj(h))
        nums.append(gexact(h, c))
        dens.append(gexact(gconj(h), c))
    B = glcm_all(dens)
    rows = [B] + [gmul(gexact(B, d), a) for a, d in zip(nums, dens)]
    g = ggcd_all(rows)
    rows = [gexact(z, g) for z in rows]
    return rows


# ------------------------------------------------------ circles and clusters
def circle_points(prime_powers, unit=(1, 0)):
    """All Gaussian integers of norm prod p^e for split p (one per unit class
    removed: returns all 4*prod(e+1) points)."""
    pts = [(1, 0)]
    for p, e in prime_powers:
        pi = split_prime_pi(p)
        pib = gconj(pi)
        new = []
        for z in pts:
            for a in range(e + 1):
                new.append(gmul(z, gmul(gpow(pi, a), gpow(pib, e - a))))
        pts = new
    out = []
    for z in pts:
        w = z
        for _ in range(4):
            out.append(w)
            w = gmul(w, (0, 1))
    return sorted(set(out))


def angle(z):
    import math
    return math.atan2(z[1], z[0]) % (2 * math.pi)


def half_cot(zi, zj, N):
    """Rational cot of half the angle from z_i to z_j (counterclockwise):
    (N + <z_i,z_j>)/det(z_i,z_j)."""
    dot = zi[0] * zj[0] + zi[1] * zj[1]
    det = zi[0] * zj[1] - zi[1] * zj[0]
    if det == 0:
        return None
    return Fraction(N + dot, det)
