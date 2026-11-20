# exact setup for q_P: primes p = 1 mod 4 in [P,2P]; fixed generators; one-unit marks (eta = 1)
from math import gcd, isqrt
from sympy import primerange
def gen(p):
    for a in range(1, isqrt(p)+1):
        b2 = p - a*a; b = isqrt(b2)
        if b*b == b2 and a % 2 == 1: return (a, b)   # a odd, b even, both >0
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def setup(P, cls=(1,4)):
    r, m = cls
    ps = [p for p in primerange(P, 2*P+1) if p % m == r]
    gs = [gen(p) for p in ps]
    q = 1
    for p in ps: q *= p
    s = len(ps)
    marks = []
    for mask in range(1 << s):          # bit k set -> conj(pi_k)
        z = (1, 0)
        for k in range(s):
            a, b = gs[k]
            z = gmul(z, (a, -b) if (mask >> k) & 1 else (a, b))
        marks.append(z)
    for z in marks: assert z[0]**2 + z[1]**2 == q
    return ps, gs, q, s, marks
def Qof(ps, mi, mj):
    d = mi ^ mj; Q = 1
    for k, p in enumerate(ps):
        if (d >> k) & 1: Q *= p
    return Q
def Tof(A, B):
    # content of (B - A)/2 ; equals |Im X| for one-unit marks
    x, y = B[0]-A[0], B[1]-A[1]
    assert x % 2 == 0 and y % 2 == 0
    return gcd(x//2, y//2)
def d2(A, B): return (B[0]-A[0])**2 + (B[1]-A[1])**2
