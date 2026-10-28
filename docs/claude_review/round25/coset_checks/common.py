"""Shared helpers for round-25 coset checks.
Normalised coordinates: x = phi/(pi/2) in R/Z.  Gap conditions for monomial v (y = exp(-w(v)/2)):
  dist(<v,x>, Z)       >= a0(y) := (2/pi) asin(y)
  dist(<v,x>, 1/2 + Z) >= a1(y) := (2/pi) asin(y/sqrt 2)
"""
import math, itertools
import numpy as np

def split_primes(upto):
    ps = []
    for n in range(5, upto + 1, 4):
        if all(n % d for d in range(2, int(n**0.5) + 1)):
            ps.append(n)
    return ps

def gaussian_prime(p):
    """a+bi with a^2+b^2=p, a>b>0 (p split prime)."""
    for b in range(1, int(math.isqrt(p)) + 1):
        a2 = p - b * b
        a = math.isqrt(a2)
        if a * a == a2 and a > b:
            return a, b
    raise ValueError(p)

def monomials(logs, T):
    """All nonzero v in Z^r with sum |v_j| logs_j <= T, one of each pair +-v (first nonzero coord > 0).
    Returns int array (n, r) and weights."""
    r = len(logs)
    out = []
    def rec(j, cur, w, started):
        if j == r:
            if started:
                out.append((tuple(cur), w))
            return
        k = 0
        while w + k * logs[j] <= T + 1e-12:
            if k == 0:
                cur.append(0); rec(j + 1, cur, w, started); cur.pop()
            else:
                for s in ((1,) if not started else (1, -1)):
                    cur.append(s * k); rec(j + 1, cur, w + k * logs[j], True); cur.pop()
            k += 1
    rec(0, [], 0.0, False)
    V = np.array([o[0] for o in out], dtype=np.int64)
    W = np.array([o[1] for o in out])
    return V, W

def halfwidths(W):
    y = np.exp(-W / 2)
    return (2 / math.pi) * np.arcsin(y), (2 / math.pi) * np.arcsin(y / math.sqrt(2))

def window_measure(W):
    a0, a1 = halfwidths(W)
    return 2 * a0 + 2 * a1   # normalised measure of both windows on R/Z

def dist_Z(t):
    return np.abs(t - np.round(t))
