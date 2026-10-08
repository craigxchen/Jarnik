"""Small exact Gaussian-integer toolkit (integers only)."""
from math import gcd, isqrt

def mul(z, w):
    return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])

def conj(z):
    return (z[0], -z[1])

def norm(z):
    return z[0]*z[0] + z[1]*z[1]

def sub(z, w):
    return (z[0]-w[0], z[1]-w[1])

def divround(z, w):
    """nearest Gaussian quotient of z/w"""
    n = norm(w)
    x = z[0]*w[0] + z[1]*w[1]
    y = z[1]*w[0] - z[0]*w[1]
    def rnd(a):  # nearest integer to a/n
        return (2*a + n) // (2*n)
    return (rnd(x), rnd(y))

def ggcd(z, w):
    while w != (0, 0):
        q = divround(z, w)
        z, w = w, sub(z, mul(q, w))
    return z

def divexact(z, w):
    n = norm(w)
    x = z[0]*w[0] + z[1]*w[1]
    y = z[1]*w[0] - z[0]*w[1]
    assert x % n == 0 and y % n == 0, (z, w)
    return (x//n, y//n)

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def rad(n):
    n = abs(n)
    r = 1
    p = 2
    while p*p <= n:
        if n % p == 0:
            r *= p
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        r *= n
    return r

def factor(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def two_squares(N):
    """all (x,y) with x^2+y^2=N (brute force over x; fine for N<=1e12 with care)"""
    pts = []
    x = 0
    r = isqrt(N)
    while x <= r:
        y2 = N - x*x
        y = isqrt(y2)
        if y*y == y2:
            for sx in ((1, -1) if x else (1,)):
                for sy in ((1, -1) if y else (1,)):
                    pts.append((sx*x, sy*y))
        x += 1
    return pts
