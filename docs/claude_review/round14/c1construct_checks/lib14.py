# Round 14 helpers (exact). Gaussian integers as tuples (x, y) = x + iy.
import math
from sympy import primerange, isprime, legendre_symbol, factorint

def chi4(q):
    return 1 if q % 4 == 1 else -1

def mq(q):
    return q - chi4(q)

def kq(q):
    """2-adic valuation of m_q (the class k(q) of round 7)."""
    m = mq(q); k = 0
    while m % 2 == 0:
        m //= 2; k += 1
    return k

def nu(p, q):
    return 1 if legendre_symbol(p % q, q) == -1 else 0

def split_primes(lo, hi):
    return [p for p in primerange(lo, hi + 1) if p % 4 == 1]

def gauss_prime(p):
    """pi = a + bi with a^2 + b^2 = p, a odd > 0, b even > 0."""
    for a in range(1, math.isqrt(p) + 1, 2):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2 and b % 2 == 0 and b > 0:
            return (a, b)
    raise ValueError(p)

def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gconj(x):
    return (x[0], -x[1])

def gpow(x, e):
    r = (1, 0)
    while e:
        if e & 1: r = gmul(r, x)
        x = gmul(x, x); e >>= 1
    return r

def Zv(zs, v):
    """prod zeta_j^{v_j} with zeta^{-1} represented by conj zeta (same class in T)."""
    z = (1, 0)
    for zz, t in zip(zs, v):
        if t:
            z = gmul(z, gpow(zz if t > 0 else gconj(zz), abs(t)))
    return z

def Ds(z, s):
    X, Y = z
    return [Y, X - Y, X, X + Y][s % 4]

# modular Gaussian arithmetic mod Q
def mmul(x, y, Q):
    return ((x[0] * y[0] - x[1] * y[1]) % Q, (x[0] * y[1] + x[1] * y[0]) % Q)

def mpow(x, e, Q):
    r = (1, 0); x = (x[0] % Q, x[1] % Q)
    while e:
        if e & 1: r = mmul(r, x, Q)
        x = mmul(x, x, Q); e >>= 1
    return r

def minv(x, Q):
    n = (x[0] * x[0] + x[1] * x[1]) % Q
    ni = pow(n, -1, Q)
    return ((x[0] * ni) % Q, (-x[1] * ni) % Q)

def tor(z, Q):
    """z / conj z mod Q (z unit mod Q)."""
    return mmul(z, minv(gconj(z), Q), Q)

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def unit_mod(s, Q):
    u = UNITS[s % 4]
    return (u[0] % Q, u[1] % Q)
