"""Shared exact helpers for the walsh-extraction checks (Gaussian integers as (a,b) int pairs)."""
import math

def is_prime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    r = int(math.isqrt(n)); f = 3
    while f <= r:
        if n % f == 0: return False
        f += 2
    return True

def split_primes(count, start=5):
    out = []; p = start
    while len(out) < count:
        if p % 4 == 1 and is_prime(p): out.append(p)
        p += 1
    return out

def gauss_prime(p):
    """Return (a,b), a>b>0, a^2+b^2=p for split prime p."""
    for b in range(1, int(math.isqrt(p)) + 1):
        a2 = p - b * b
        a = math.isqrt(a2)
        if a * a == a2 and a > 0:
            return (max(a, b), min(a, b))
    raise ValueError(p)

def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gconj(x):
    return (x[0], -x[1])

def gpow(x, e):
    r = (1, 0)
    for _ in range(e): r = gmul(r, x)
    return r

def gnorm(x):
    return x[0] * x[0] + x[1] * x[1]

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def ggcd(x, y):
    while y != (0, 0):
        # division with rounding in Q(i)
        n = gnorm(y)
        num = gmul(x, gconj(y))
        q = ((2 * num[0] + n) // (2 * n), (2 * num[1] + n) // (2 * n))
        r = (x[0] - gmul(q, y)[0], x[1] - gmul(q, y)[1])
        x, y = y, r
    return x

def gdivexact(x, y):
    n = gnorm(y); num = gmul(x, gconj(y))
    assert num[0] % n == 0 and num[1] % n == 0, (x, y)
    return (num[0] // n, num[1] // n)
