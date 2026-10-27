"""Common regime parameters for round-16 antitheory checks (T = 1).
d = n-1, w = sqrt(8 d ln 4r), j = floor(w/2), m' = ceil(4w),
m_max = largest m with 8 binom(m+d-1, d) <= 3^d (exact big integers when feasible),
W_min = 3 n ln n - 14 n (two_thirds (T3), C <= 1), abar = W/(2d), eps_0 = d e^{-W/4}/(m'+1) (C = 1).
"""
import math
from mpmath import mp, mpf, log, sqrt, ceil, floor, exp
mp.dps = 50

def wjm(n, r):
    d = n - 1
    w = sqrt(8 * d * log(4 * mpf(r)))
    j = int(floor(w / 2))
    m = int(ceil(4 * w))
    return d, w, j, m

def mmax_ok(n, m):
    """exact test 8 binom(m+d-1,d) <= 3^d (R1) for moderate n; returns None if too big."""
    d = n - 1
    if d > 300000:
        # log test with mpmath (50 digits) -- evidence-level only for huge n
        lhs = log(8) + mp.loggamma(m + d) - mp.loggamma(d + 1) - mp.loggamma(m)
        return bool(lhs <= d * log(3))
    return 8 * math.comb(m + d - 1, d) <= 3 ** d

def Wmin(n):
    return 3 * n * log(n) - 14 * n

def ln_inv_eps0(n, m, W):
    d = n - 1
    return W / 4 - log(mpf(d) / (m + 1))
