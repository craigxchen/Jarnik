"""Exact Gaussian-integer helpers (pure integers, no floating point in identities)."""
from math import gcd, isqrt
import random

def mul(z, w):
    return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])

def add(z, w):
    return (z[0]+w[0], z[1]+w[1])

def sub(z, w):
    return (z[0]-w[0], z[1]-w[1])

def conj(z):
    return (z[0], -z[1])

def norm(z):
    return z[0]*z[0] + z[1]*z[1]

def smul(c, z):
    return (c*z[0], c*z[1])

def gpow(z, e):
    r = (1, 0)
    for _ in range(e):
        r = mul(r, z)
    return r

def divround(z, w):
    n = norm(w)
    x = z[0]*w[0] + z[1]*w[1]
    y = z[1]*w[0] - z[0]*w[1]
    def rnd(a):
        return (2*a + n) // (2*n)
    return (rnd(x), rnd(y))

def ggcd(z, w):
    while w != (0, 0):
        q = divround(z, w)
        z, w = w, sub(z, mul(q, w))
    return z

def divides(w, z):
    n = norm(w)
    x = z[0]*w[0] + z[1]*w[1]
    y = z[1]*w[0] - z[0]*w[1]
    return x % n == 0 and y % n == 0

def divexact(z, w):
    n = norm(w)
    x = z[0]*w[0] + z[1]*w[1]
    y = z[1]*w[0] - z[0]*w[1]
    assert x % n == 0 and y % n == 0, (z, w)
    return (x//n, y//n)

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def content(z):
    return gcd(abs(z[0]), abs(z[1]))

def primitive_direction(z):
    c = content(z)
    return (z[0]//c, z[1]//c)

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f = 3
    while f*f <= n:
        if n % f == 0:
            return False
        f += 2
    return True

def split_primes(limit):
    return [p for p in range(5, limit) if p % 4 == 1 and is_prime(p)]

def gaussian_prime_above(p):
    for a in range(1, isqrt(p)+1):
        b2 = p - a*a
        b = isqrt(b2)
        if b*b == b2 and b > 0:
            return (a, b)
    raise ValueError(p)

def vp(n, p):
    n = abs(n)
    if n == 0:
        return 10**9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v
