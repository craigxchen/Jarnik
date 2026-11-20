"""Exact Gaussian-integer utilities (referee's own code)."""
from math import isqrt, gcd
from sympy import factorint, primerange

def mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])

def conj(z):
    return (z[0], -z[1])

def nrm(z):
    return z[0]*z[0]+z[1]*z[1]

def powg(z, k):
    r = (1, 0)
    for _ in range(k):
        r = mul(r, z)
    return r

def divexact(z, w):
    """z/w if exact in Z[i], else None"""
    n = nrm(w)
    a, b = mul(z, conj(w))
    if a % n or b % n:
        return None
    return (a//n, b//n)

def ggcd(z, w):
    while w != (0, 0):
        n = nrm(w)
        a, b = mul(z, conj(w))
        # rounded quotient
        qa = (2*a + n)//(2*n); qb = (2*b + n)//(2*n)
        q = (qa, qb)
        qw = mul(q, w)
        z, w = w, (z[0]-qw[0], z[1]-qw[1])
    return z

_pi_cache = {}
def gprime(p):
    """a Gaussian prime pi with pi*conj(pi)=p, p = 1 mod 4"""
    if p in _pi_cache:
        return _pi_cache[p]
    for a in range(1, isqrt(p)+1):
        b2 = p - a*a
        b = isqrt(b2)
        if b*b == b2:
            _pi_cache[p] = (a, b)
            return (a, b)
    raise ValueError(p)

def val(z, pi):
    k = 0
    while True:
        q = divexact(z, pi)
        if q is None:
            return k
        z = q; k += 1

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def nu_formula(p, D):
    """number of x mod p with x^2+y^2 = D solvable mod p, p odd prime"""
    D %= p
    chi4 = 1 if p % 4 == 1 else -1
    if D == 0:
        return p if p % 4 == 1 else 1
    leg = pow(D, (p-1)//2, p)
    leg = 1 if leg == 1 else -1
    return (p - chi4 + 1 + leg)//2

def nu_brute(p, D):
    sq = {x*x % p for x in range(p)}
    return sum(((D - x*x) % p) in sq for x in range(p))
