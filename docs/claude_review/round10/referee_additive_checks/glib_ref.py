"""Independent Gaussian-integer toolkit for the referee checks (not derived from additive_checks/glib.py).

Gaussian integers are pairs (a, b) = a + b i of Python ints. All arithmetic exact.
"""
import math, random
from sympy import isprime, factorint
import mpmath as mp

mp.mp.dps = 60


def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gconj(x):
    return (x[0], -x[1])


def gnorm(x):
    return x[0] * x[0] + x[1] * x[1]


def gsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def gdivmod(x, y):
    """Euclidean division x = q y + r with Norm r < Norm y (rounded quotient)."""
    n = gnorm(y)
    num = gmul(x, gconj(y))
    # round num/n to nearest Gaussian integer
    q = ((2 * num[0] + n) // (2 * n), (2 * num[1] + n) // (2 * n))
    r = gsub(x, gmul(q, y))
    return q, r


def gexact_div(x, y):
    n = gnorm(y)
    num = gmul(x, gconj(y))
    if num[0] % n or num[1] % n:
        return None
    return (num[0] // n, num[1] // n)


def ggcd(x, y):
    while y != (0, 0):
        _, r = gdivmod(x, y)
        x, y = y, r
    return x


def gprime_above(p):
    """A Gaussian prime pi with Norm pi = p, p = 1 mod 4."""
    assert p % 4 == 1
    for g in range(2, p):
        t = pow(g, (p - 1) // 4, p)
        if t * t % p == p - 1:
            break
    pi = ggcd((p, 0), (t, 1))
    assert gnorm(pi) == p
    return pi


def gval(z, pi):
    """pi-adic valuation of a nonzero Gaussian integer z."""
    v = 0
    while True:
        q = gexact_div(z, pi)
        if q is None:
            return v
        z = q
        v += 1


UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]


def points_on_circle(fac):
    """All lattice points with Norm = N = prod p^e, all p = 1 mod 4 (dict fac)."""
    pts = [(1, 0)]
    for p, e in fac.items():
        pi = gprime_above(p)
        pib = gconj(pi)
        new = []
        pw = [(1, 0)]
        for _ in range(e):
            pw.append(gmul(pw[-1], pi))
        pwb = [(1, 0)]
        for _ in range(e):
            pwb.append(gmul(pwb[-1], pib))
        facs = [gmul(pw[a], pwb[e - a]) for a in range(e + 1)]
        for z in pts:
            for f in facs:
                new.append(gmul(z, f))
        pts = new
    out = set()
    for z in pts:
        for u in UNITS:
            out.add(gmul(z, u))
    return sorted(out)


def angle(z):
    return mp.atan2(z[1], z[0]) % (2 * mp.pi)


def windows(pts_sorted_with_angles, D):
    """Maximal windows of consecutive points (circularly) of angular width <= D. Returns index lists."""
    P = pts_sorted_with_angles
    n = len(P)
    res = []
    prev_end = -1
    j = 0
    for i in range(n):
        if j < i:
            j = i
        while j + 1 < i + n:
            dth = P[(j + 1) % n][1] - P[i][1]
            if (j + 1) >= n:
                dth += 2 * mp.pi
            if dth <= D:
                j += 1
            else:
                break
        if j > prev_end:
            res.append([(k % n) for k in range(i, j + 1)])
            prev_end = j
    return res


def cluster_gcd(zs):
    g = zs[0]
    for z in zs[1:]:
        g = ggcd(g, z)
    return g


def normalise(zs):
    g = cluster_gcd(zs)
    out = [gexact_div(z, g) for z in zs]
    return out, g


def profile(zs, N):
    """Exponent profile of a normalised cluster on Norm = N. Returns (primes, e, pis, a[x][j])."""
    fac = factorint(N)
    primes = sorted(fac)
    pis = [gprime_above(p) for p in primes]
    e = [fac[p] for p in primes]
    a = [[gval(z, pi) for pi in pis] for z in zs]
    return primes, e, pis, a


def random_split_N(rng, maxN, minN=10**4):
    """Random odd N with all prime factors = 1 mod 4."""
    smallp = [p for p in range(5, 2000) if p % 4 == 1 and isprime(p)]
    while True:
        N = 1
        fac = {}
        k = rng.randint(2, 7)
        for _ in range(k):
            p = rng.choice(smallp[:rng.choice([8, 20, 60, len(smallp)])])
            N *= p
            fac[p] = fac.get(p, 0) + 1
            if N > maxN:
                break
        if minN <= N <= maxN:
            return N, fac
