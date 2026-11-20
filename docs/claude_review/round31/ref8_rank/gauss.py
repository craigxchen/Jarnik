"""Shared exact helpers (own code, referee of Part VIII.1).

Gaussian integers are pairs (x, y) of Python ints.  All rank computations are exact
(fraction-free elimination over Z).
"""
from math import atan2, log, sqrt, pi, isqrt
from itertools import product
from sympy import factorint, primerange


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def gpow(z, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, z)
    return r


_gp_cache = {}


def gauss_prime(p):
    """a+bi with a^2+b^2=p, a>b>0 (p = 1 mod 4)."""
    if p in _gp_cache:
        return _gp_cache[p]
    assert p % 4 == 1
    # Cornacchia / Hermite-Serret via sqrt(-1) mod p
    for g in range(2, p):
        t = pow(g, (p - 1) // 4, p)
        if t * t % p == p - 1:
            break
    a, b = p, t
    lim = isqrt(p)
    while b > lim:
        a, b = b, a % b
    c = isqrt(p - b * b)
    assert b * b + c * c == p
    x, y = max(b, c), min(b, c)
    _gp_cache[p] = (x, y)
    return (x, y)


def rank_int(rows):
    """Exact rank of an integer matrix (list of lists)."""
    m = [list(r) for r in rows]
    if not m:
        return 0
    ncol = len(m[0])
    rk = 0
    prow = 0
    for c in range(ncol):
        piv = None
        for r in range(prow, len(m)):
            if m[r][c] != 0:
                piv = r
                break
        if piv is None:
            continue
        m[prow], m[piv] = m[piv], m[prow]
        pv = m[prow][c]
        for r in range(prow + 1, len(m)):
            if m[r][c] != 0:
                f = m[r][c]
                m[r] = [pv * m[r][k] - f * m[prow][k] for k in range(ncol)]
                # keep entries small
                from math import gcd
                g = 0
                for v in m[r]:
                    g = gcd(g, v)
                if g > 1:
                    m[r] = [v // g for v in m[r]]
        prow += 1
        rk += 1
        if prow == len(m):
            break
    return rk


def affine_rank(rows):
    if len(rows) <= 1:
        return 0
    r0 = rows[0]
    return rank_int([[a - b for a, b in zip(r, r0)] for r in rows[1:]])


PRIMES3 = [l for l in primerange(3, 200000) if l % 4 == 3]
_cum = []
_s = 0.0
for _l in PRIMES3:
    _s += log(_l) / (_l + 1)
    _cum.append(_s)


def initial_S(logC, extra=0.0):
    """Shortest initial segment S of primes 3 mod 4 with sum log l/(l+1) + extra > logC."""
    if extra > logC:
        return []
    for k, s in enumerate(_cum):
        if s + extra > logC:
            return PRIMES3[:k + 1]
    raise ValueError("C too large")


def indicator_block(points, mod_fn):
    """One-hot indicator of residue classes; returns list of rows and number of classes."""
    keys = [mod_fn(z) for z in points]
    cls = sorted(set(keys))
    idx = {k: i for i, k in enumerate(cls)}
    rows = []
    for k in keys:
        v = [0] * len(cls)
        v[idx[k]] = 1
        rows.append(v)
    return rows, len(cls)
