"""Exact Gaussian-integer helpers (pairs of Python ints). Round 27, lane stab."""
import random
from math import isqrt
from itertools import combinations

def gmul(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])

def gconj(x):
    return (x[0], -x[1])

def gnorm(x):
    return x[0]*x[0] + x[1]*x[1]

def gsub(x, y):
    return (x[0]-y[0], x[1]-y[1])

def gprod(lst):
    r = (1, 0)
    for z in lst:
        r = gmul(r, z)
    return r

def gpow(x, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, x)
    return r

def gdivmod(x, y):
    # x / y rounded
    n = gnorm(y)
    num = gmul(x, gconj(y))
    def rnd(a):
        # round a/n to nearest integer
        return (2*a + n) // (2*n)
    q = (rnd(num[0]), rnd(num[1]))
    r = gsub(x, gmul(q, y))
    return q, r

def gexactdiv(x, y):
    n = gnorm(y)
    num = gmul(x, gconj(y))
    assert num[0] % n == 0 and num[1] % n == 0, (x, y)
    return (num[0]//n, num[1]//n)

def gdivides(y, x):
    n = gnorm(y)
    num = gmul(x, gconj(y))
    return num[0] % n == 0 and num[1] % n == 0

def ggcd(x, y):
    while y != (0, 0):
        _, r = gdivmod(x, y)
        x, y = y, r
    return x

def ggcd_list(lst):
    g = (0, 0)
    for z in lst:
        g = ggcd(g, z)
    return g

def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n-1):
            continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1:
                break
        else:
            return False
    return True

def gauss_prime(p):
    """pi = a+bi with a^2+b^2=p, a>b>0 (fixed choice)."""
    assert p % 4 == 1
    # find sqrt(-1) mod p
    for c in range(2, p):
        if pow(c, (p-1)//2, p) == p-1:
            r = pow(c, (p-1)//4, p)
            break
    a, b = p, r
    lim = isqrt(p)
    while b > lim:
        a, b = b, a % b
    c = isqrt(p - b*b)
    assert b*b + c*c == p
    x, y = max(b, c), min(b, c)
    return (x, y)

def random_split_primes(r, lo, hi, rng):
    out = set()
    while len(out) < r:
        p = rng.randrange(lo, hi)
        if p % 4 == 1 and is_prime(p):
            out.add(p)
    return sorted(out)

def uniquely_divides_count(z, pi):
    k = 0
    while gdivides(pi, z):
        z = gexactdiv(z, pi); k += 1
    return k
