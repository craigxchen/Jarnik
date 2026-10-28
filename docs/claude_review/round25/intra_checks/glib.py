"""Exact Gaussian-integer toolkit for round 25 intra-cluster checks.
Gaussian integers are tuples (a,b); Gaussian rationals are tuples (Fraction, Fraction)."""
import math
from fractions import Fraction as Fr
from math import isqrt
from sympy import factorint

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod(a, b):
    n = gnorm(b); c = gmul(a, gconj(b))
    q = ((2*c[0]+n)//(2*n), (2*c[1]+n)//(2*n))
    qb = gmul(q, b); return q, (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod(a, b); a, b = b, r
    return a
def gexdiv(a, b):
    q, r = gdivmod(a, b); assert r == (0, 0), (a, b); return q
def gpow(a, e):
    r = (1, 0)
    for _ in range(e): r = gmul(r, a)
    return r

# Gaussian rationals
def qmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def qinv(a):
    n = a[0]*a[0]+a[1]*a[1]; return (a[0]/n, -a[1]/n)
def qadd(a, b): return (a[0]+b[0], a[1]+b[1])
def qnorm(a): return a[0]*a[0]+a[1]*a[1]
def q_of(g): return (Fr(g[0]), Fr(g[1]))

def points(N):
    P = []
    for x in range(-isqrt(N), isqrt(N)+1):
        y2 = N-x*x; y = isqrt(y2)
        if y*y == y2:
            P.append((x, y))
            if y: P.append((x, -y))
    P.sort(key=lambda p: math.atan2(p[1], p[0])); return P

def best_window(N, M):
    P = points(N); n = len(P); best = None
    for i in range(n):
        W = [P[(i+j) % n] for j in range(M)]
        d = max(math.dist(a, b) for a in W for b in W)
        if best is None or d < best[0]: best = (d, W)
    return best

def split_prime_pi(p):
    for a in range(1, isqrt(p)+1):
        b2 = p-a*a; b = isqrt(b2)
        if b*b == b2: return (a, b)
    raise ValueError(p)

def vpi(z, pi):
    """valuation of Gaussian integer z at Gaussian prime pi (z != 0)."""
    v = 0
    while True:
        q, r = gdivmod(z, pi)
        if r != (0, 0): return v
        z = q; v += 1

def normalised_cluster(N, M):
    d, W = best_window(N, M)
    g = W[0]
    for z in W[1:]: g = ggcd(g, z)
    Z = [gexdiv(z, g) for z in W]
    Nn = N // gnorm(g)
    return d, Z, Nn, g

def profile(Z, Nn):
    """returns list of (p, pi, a_vector) for split primes p | Nn; a_x = v_pi(z_x)."""
    out = []
    for p, e in factorint(Nn).items():
        assert p % 4 == 1, p
        pi = split_prime_pi(p)
        out.append((p, e, pi, [vpi(z, pi) for z in Z]))
    return out

def char_value(Z, k):
    """Z^k for k in Z^M with sum 0, exact Gaussian rational."""
    assert sum(k) == 0
    num = (1, 0); den = (1, 0)
    for z, kk in zip(Z, k):
        if kk > 0: num = gmul(num, gpow(z, kk))
        elif kk < 0: den = gmul(den, gpow(z, -kk))
    return qmul(q_of(num), qinv(q_of(den)))

def width(supp, v):
    vals = [sum(a*b for a, b in zip(k, v)) for k in supp]
    return max(vals)-min(vals)
