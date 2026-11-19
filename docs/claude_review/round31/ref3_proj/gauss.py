# Independent Gaussian-integer and star-matrix tools for the ref3_proj report.
from math import gcd, isqrt
from fractions import Fraction
import itertools, sympy

def gm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gc(a): return (a[0], -a[1])
def gn(a): return a[0]*a[0]+a[1]*a[1]
def gdiv_exact(a, b):
    n = gn(b); num = gm(a, gc(b))
    if num[0] % n or num[1] % n: return None
    return (num[0]//n, num[1]//n)
def gmod(a, b):
    n = gn(b); num = gm(a, gc(b))
    q = ((2*num[0]+n)//(2*n), (2*num[1]+n)//(2*n))
    qb = gm(q, b); return (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0):
        a, b = b, gmod(a, b)
    return a
def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gm(r, a)
    return r
def split_prime(p):
    r = sympy.sqrt_mod(-1, p)
    pi = ggcd((p, 0), (r, 1))
    assert gn(pi) == p
    return pi
def val(z, pi):
    v = 0
    while True:
        q = gdiv_exact(z, pi)
        if q is None: return v
        z = q; v += 1
def det(a, b): return a[0]*b[1]-a[1]*b[0]   # [a,b] = Im(conj(a) b)

UNITW = {(1, 0): (1, 0), (0, 1): (1, 1), (-1, 0): (0, 1), (0, -1): (0, -1)}
# w with w/conj(w) = unit: 1->1, i->1+i, -1->i, -i->-i*(1+i)/... check below

def unit_w(u):
    for w in [(1, 0), (1, 1), (0, 1), (1, -1)]:
        q = gdiv_exact(w, gc(w))
        if q == u: return w
    raise ValueError(u)

def stars(points):
    """points: list of Gaussian integers of a common norm. Returns dict with
    primes, exponents e_p, arrays a[i][p], stars P_i' (unit-absorbed), g, d, t."""
    g0 = (0, 0)
    for z in points: g0 = ggcd(g0, z)
    zs = [gdiv_exact(z, g0) for z in points]
    N = gn(zs[0]); assert all(gn(z) == N for z in zs)
    fac = sympy.factorint(N)
    primes = sorted(fac)
    for p in primes: assert p % 4 == 1, ("non-split prime left", p)
    pis = {p: split_prime(p) for p in primes}
    m = len(zs)
    a = [[val(z, pis[p]) for p in primes] for z in zs]
    units = []
    Ps = []
    for i, z in enumerate(zs):
        base = (1, 0); P = (1, 0)
        for k, p in enumerate(primes):
            e = fac[p]; ai = a[i][k]
            assert val(z, gc(pis[p])) == e - ai
            base = gm(base, gm(gpow(pis[p], ai), gpow(gc(pis[p]), e-ai)))
            P = gm(P, gpow(pis[p], ai))
        u = gdiv_exact(z, base); assert gn(u) == 1
        units.append(u); Ps.append(P)
    # absorb units relative to first point
    Pp = []
    for i in range(m):
        rel = gdiv_exact(units[i], units[0])
        w = unit_w(rel)
        Pp.append(gm(w, Ps[i]))
    g = {}; d = {}; t = {}
    for i in range(m):
        for j in range(i+1, m):
            gg = 1; dd = 1
            for k, p in enumerate(primes):
                gg *= p**min(a[i][k], a[j][k]); dd *= p**abs(a[i][k]-a[j][k])
            D = det(Pp[i], Pp[j])
            assert D % gg == 0, "g does not divide minor"
            g[i, j] = gg; d[i, j] = dd; t[i, j] = D//gg
            assert t[i, j] != 0
    return dict(zs=zs, N=N, fac=fac, primes=primes, a=a, P=Pp, g=g, d=d, t=t, m=m)
