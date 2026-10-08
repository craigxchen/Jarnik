"""Referee's independent helpers (no code shared with the author's gi.py)."""
import math, itertools, random

def gmul(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])

def gconj(a):
    return (a[0], -a[1])

def gpow(a, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, a)
    return r

def gsub(a, b):
    return (a[0]-b[0], a[1]-b[1])

def gadd(a, b):
    return (a[0]+b[0], a[1]+b[1])

def gscal(c, a):
    return (c*a[0], c*a[1])

def gnorm(a):
    return a[0]*a[0] + a[1]*a[1]

def gdivexact(a, b):
    """a/b in Z[i], asserting exactness."""
    n = gnorm(b)
    num = gmul(a, gconj(b))
    assert num[0] % n == 0 and num[1] % n == 0, (a, b)
    return (num[0]//n, num[1]//n)

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def two_squares(p):
    """p = 1 mod 4 prime -> (a,b) with a^2+b^2=p, a odd > 0, b even > 0 (brute force, p small)."""
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a*a
        b = math.isqrt(b2)
        if b*b == b2 and b > 0:
            if a % 2 == 1:
                return (a, b)
            return (b, a)
    raise ValueError(p)

def is_prime(n):
    if n < 2:
        return False
    for q in range(2, int(math.isqrt(n)) + 1):
        if n % q == 0:
            return False
    return True

SPLIT = [p for p in range(5, 400) if p % 4 == 1 and is_prime(p)]

def factor(n):
    n = abs(n)
    out = {}
    q = 2
    while q*q <= n:
        while n % q == 0:
            out[q] = out.get(q, 0) + 1
            n //= q
        q += 1 if q == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def vp(n, p):
    n = abs(n)
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def det2(u, v):
    return u[0]*v[1] - u[1]*v[0]
