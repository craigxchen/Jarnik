"""Independent Gaussian-integer toolkit for the kill report (does not import the proposal's code)."""
import math
from math import isqrt
from sympy import factorint, sqrt_mod

def gmul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])

def gsub(a, b):
    return (a[0]-b[0], a[1]-b[1])

def gconj(a):
    return (a[0], -a[1])

def gnorm(a):
    return a[0]*a[0]+a[1]*a[1]

def gdivmod_round(a, b):
    # nearest-integer quotient
    n = gnorm(b)
    num = gmul(a, gconj(b))
    q = ((2*num[0] + n)//(2*n), (2*num[1] + n)//(2*n))
    r = gsub(a, gmul(q, b))
    return q, r

def gdiv_exact(a, b):
    n = gnorm(b)
    num = gmul(a, gconj(b))
    assert num[0] % n == 0 and num[1] % n == 0, (a, b)
    return (num[0]//n, num[1]//n)

def gdivides(b, a):
    n = gnorm(b)
    num = gmul(a, gconj(b))
    return num[0] % n == 0 and num[1] % n == 0

def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod_round(a, b)
        a, b = b, r
    return a

def gprime_above(p):
    """Gaussian prime pi = x+iy with x^2+y^2 = p, p = 1 mod 4."""
    if p == 2:
        return (1, 1)
    r = sqrt_mod(-1, p)
    a, b = (p, 0), (r, 1)
    g = ggcd(a, b)
    assert gnorm(g) == p, (p, g)
    return g

def gval(z, pi):
    """valuation of Gaussian integer z at Gaussian prime pi (z != 0)."""
    assert z != (0, 0)
    v = 0
    while gdivides(pi, z):
        z = gdiv_exact(z, pi)
        v += 1
    return v

def lattice_points(n):
    out = []
    for x in range(-isqrt(n), isqrt(n)+1):
        y2 = n - x*x
        y = isqrt(y2)
        if y*y == y2:
            out.append((x, y))
            if y:
                out.append((x, -y))
    return out

def ang(z):
    return math.atan2(z[1], z[0])
