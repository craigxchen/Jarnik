# Independent referee library (does not import the author's glib.py).
# Gaussian integers as Python complex-free tuples (a,b); exact Q(i) via Fraction pairs.
from fractions import Fraction as Fr
from math import isqrt, log, atan2, pi, gcd
from itertools import combinations
import random

def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gconj(x): return (x[0], -x[1])
def gnorm(x): return x[0]*x[0]+x[1]*x[1]
def gdivmod(x, y):
    n = gnorm(y); num = gmul(x, gconj(y))
    q = ((2*num[0]+n)//(2*n), (2*num[1]+n)//(2*n))
    r = (x[0]-(q[0]*y[0]-q[1]*y[1]), x[1]-(q[0]*y[1]+q[1]*y[0]))
    return q, r
def gdiv_exact(x, y):
    q, r = gdivmod(x, y); assert r == (0, 0); return q
def ggcd(x, y):
    while y != (0, 0):
        _, r = gdivmod(x, y); x, y = y, r
    return x
def divides(d, x):
    _, r = gdivmod(x, d); return r == (0, 0)

def factor(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def split_prime(p):
    # pi with N(pi)=p, p = 1 mod 4
    for a in range(1, isqrt(p)+1):
        b2 = p-a*a; b = isqrt(b2)
        if b*b == b2: return (a, b)
    raise ValueError

def lattice_points(N):
    pts = []
    for x in range(-isqrt(N), isqrt(N)+1):
        y2 = N-x*x; y = isqrt(y2)
        if y*y == y2:
            pts.append((x, y))
            if y: pts.append((x, -y))
    return sorted(set(pts), key=lambda z: atan2(z[1], z[0]) % (2*pi))

def vpi(z, p_pi):
    v = 0
    while divides(p_pi, z): z = gdiv_exact(z, p_pi); v += 1
    return v

# exact Q(i) arithmetic: element = (Fr re, Fr im)
def qmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qinv(x):
    n = x[0]*x[0]+x[1]*x[1]; return (x[0]/n, -x[1]/n)
def qadd(x, y): return (x[0]+y[0], x[1]+y[1])
def qpow(x, k):
    if k < 0: x = qinv(x); k = -k
    r = (Fr(1), Fr(0))
    for _ in range(k): r = qmul(r, x)
    return r
def tofr(z): return (Fr(z[0]), Fr(z[1]))
def qnorm(x): return x[0]*x[0]+x[1]*x[1]

def char_val(pts, k):
    # prod pts[j]^k[j] exactly in Q(i)
    r = (Fr(1), Fr(0))
    for j, e in enumerate(k):
        if e: r = qmul(r, qpow(tofr(pts[j]), e))
    return r

# exact rank over Q (fraction-free Gaussian elimination via Fractions)
def rank_Q(rows):
    A = [[Fr(v) for v in r] for r in rows]
    if not A: return 0
    m, n = len(A), len(A[0]); rk = 0; col = 0
    for col in range(n):
        piv = None
        for i in range(rk, m):
            if A[i][col] != 0: piv = i; break
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(m):
            if i != rk and A[i][col] != 0:
                f = A[i][col]/A[rk][col]
                A[i] = [a-f*b for a, b in zip(A[i], A[rk])]
        rk += 1
        if rk == m: break
    return rk

def monos(dim, deg):
    # exponent vectors of total degree < deg in dim variables
    out = []
    def rec(i, left, cur):
        if i == dim: out.append(tuple(cur)); return
        for e in range(left+1): rec(i+1, left-e, cur+[e])
    for d in range(deg):
        rec(0, d, []) if False else None
    # simpler: all with sum < deg
    res = []
    def rec2(i, left, cur):
        if i == dim: res.append(tuple(cur)); return
        for e in range(left+1): rec2(i+1, left-e, cur+[e])
    rec2(0, deg-1, []) if deg >= 1 else None
    return res

def reg(S):
    # S: list of k in Z^M with sum 0. reg = max m such that some nonzero c on S has sum c_k P(k)=0 for all deg P<m,
    # with P polynomial in coordinates k_1..k_{M-1} (sum-zero slice).  reg(S) = max{m : rank Poly_{<m}|_S < |S|}.
    S = [tuple(k) for k in S]; M = len(S[0]); d = M-1
    m = 0
    while True:
        mons = monos(d, m+1)  # degree < m+1
        rows = [[_mon(k, e) for k in S] for e in mons]
        if rank_Q(rows) == len(S): return m
        m += 1
def _mon(k, e):
    r = 1
    for ki, ei in zip(k, e): r *= ki**ei
    return r

def width(S, v):
    vals = [sum(a*b for a, b in zip(k, v)) for k in S]
    return max(vals)-min(vals)
