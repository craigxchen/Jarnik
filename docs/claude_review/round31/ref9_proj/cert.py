"""Normalised cluster, exponent rows, certificates T(c), b, projections. Referee's own code."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9_proj')
from gz import *
from math import atan2, pi, sqrt, log, sin, cos
from fractions import Fraction
from itertools import combinations, product
from sympy import factorint, primerange

ETAS = [(1, 0), (0, 1), (1, 1), (1, -1)]

def normalise(pts):
    g = (0, 0)
    for z in pts:
        g = ggcd(g, z) if g != (0, 0) else z
    zs = [divexact(z, g) for z in pts]
    Np = nrm(zs[0])
    assert all(nrm(z) == Np for z in zs)
    return g, zs, Np

def rows_of(zs, Np):
    f = factorint(Np)
    primes = sorted(f)
    assert all(p % 4 == 1 for p in primes), (Np, f)
    pis = [gprime(p) for p in primes]
    rows = [[val(z, pi_) for pi_ in pis] for z in zs]
    e = [f[p] for p in primes]
    for r, z in zip(rows, zs):
        for a, ep, pi_ in zip(r, e, pis):
            assert val(z, conj(pi_)) == ep - a
    return primes, pis, e, rows

def certificate(zs, Np, primes, pis, e, rows, c):
    """returns (T, b, u) with prod z^c / N' = u alpha/conj(alpha), b = eta alpha, b/bbar = prod z^c/N'.
    Verified exactly."""
    assert sum(c) == 2
    v = [sum(ci*r[k] for ci, r in zip(c, rows)) - e[k] for k in range(len(primes))]
    T = 1; alpha = (1, 0)
    for p, pi_, vk in zip(primes, pis, v):
        T *= p**abs(vk)
        alpha = mul(alpha, powg(pi_ if vk > 0 else conj(pi_), abs(vk)))
    assert nrm(alpha) == T
    # numerator: prod_{c>0} z^c prod_{c<0} conj(z)^{|c|};  prod z^c = num / N'^{neg}
    num = (1, 0); neg = 0
    for ci, z in zip(c, zs):
        if ci > 0: num = mul(num, powg(z, ci))
        elif ci < 0: num = mul(num, powg(conj(z), -ci)); neg += -ci
    # want num * conj(alpha) == u * alpha * N'^{1+neg}
    lhs = mul(num, conj(alpha))
    Nk = Np**(1+neg)
    us = [u for u in UNITS if mul(mul(u, alpha), (Nk, 0)) == lhs]
    assert len(us) == 1, "factorisation identity failed"
    u = us[0]
    etas = [eta for eta in ETAS if mul(eta, (1, 0)) == mul(u, conj(eta))]
    assert len(etas) == 1
    b = mul(etas[0], alpha)
    # check b/bbar = prod z^c / N'  <=>  num * conj(b) == b * N'^{1+neg}
    assert mul(num, conj(b)) == mul(b, (Nk, 0))
    return T, b, u, v

def projections(zs, b):
    return sorted({mul(conj(b), z)[0] for z in zs})

def angles(zs):
    a = [atan2(z[1], z[0]) for z in zs]
    # unwrap around first
    a0 = a[0]
    a = [a0 + ((t - a0 + pi) % (2*pi) - pi) for t in a]
    lo, hi = min(a), max(a)
    return a, (lo+hi)/2, hi-lo

def sieve_check(xs, D):
    """exact larger sieve inequality: sum_{p<=M} log p (n^2/nu_p - n) <= sum_{i!=j} log|x_i-x_j| <= n(n-1) log X"""
    n = len(xs)
    if n < 2: return None
    X = xs[-1]-xs[0]
    M = 2*n  # we will call with M passed separately
    return X

def sieve_terms(xs, D, M):
    n = len(xs)
    X = xs[-1]-xs[0]
    lhs = 0.0
    for p in primerange(2, M+1):
        if p == 2:
            nu = len({x % 2 for x in range(2) for y in range(2) if (x*x+y*y-D) % 2 == 0})
        else:
            nu = nu_formula(p, D)
        lhs += log(p)*(n*n/nu - n)
    mid = sum(log(abs(a-b)) for a in xs for b in xs if a != b)
    rhs = n*(n-1)*log(X) if X > 0 else 0.0
    return lhs, mid, rhs

def bad_sum(Np, M):
    return sum(log(p)/p for p in factorint(Np) if p <= M and p % 4 == 1)
