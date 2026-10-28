# Referee J4: exact Gaussian-integer helpers (Python ints only).
from math import isqrt, atan2
from fractions import Fraction as Fr

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def norm(a): return a[0]*a[0]+a[1]*a[1]
def gdivexact(a, b):
    """a/b if b | a in Z[i], else None."""
    nb = norm(b); c = gmul(a, gconj(b))
    if c[0] % nb or c[1] % nb: return None
    return (c[0]//nb, c[1]//nb)
def gmod(a, b):
    nb = norm(b); c = gmul(a, gconj(b))
    # nearest-integer rounding with exact rationals
    q = (round(Fr(c[0], nb)), round(Fr(c[1], nb)))
    qb = gmul(q, b)
    return (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0): a, b = b, gmod(a, b)
    return a
def ggcd_list(L):
    g = (0, 0)
    for z in L: g = ggcd(g, z)
    return g
def val(z, pi):
    c = 0
    while True:
        q = gdivexact(z, pi)
        if q is None: return c
        z = q; c += 1
def factor_int(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f
def gauss_prime_above(p):
    """pi = a+bi with a^2+b^2=p, p = 1 mod 4 (fixed choice a>b>0)."""
    for a in range(1, isqrt(p)+1):
        b2 = p-a*a; b = isqrt(b2)
        if b*b == b2 and a > b > 0: return (a, b)
    raise ValueError(p)
def lattice_points(N):
    pts = set()
    for x in range(0, isqrt(N)+1):
        y2 = N-x*x; y = isqrt(y2)
        if y*y == y2:
            for sx in (1, -1):
                for sy in (1, -1): pts.add((sx*x, sy*y))
    return sorted(pts, key=lambda p: atan2(p[1], p[0]))
def is_unit(a): return norm(a) == 1
