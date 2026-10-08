# aux_gauss.py -- exact Gaussian-integer toolkit for the auxiliary-methods write-up (auxiliary.md).
# All arithmetic on Gaussian integers is exact (Python ints).  Floats are used ONLY to sort points
# by angle and to locate candidate clusters; every claimed identity is then re-checked exactly.
import math, random
from itertools import combinations
from fractions import Fraction

# ---------- Gaussian integers as (a,b) = a+bi ----------
def gadd(x, y): return (x[0] + y[0], x[1] + y[1])
def gsub(x, y): return (x[0] - y[0], x[1] - y[1])
def gmul(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def gconj(x): return (x[0], -x[1])
def gnorm(x): return x[0] * x[0] + x[1] * x[1]
def gpow(x, n):
    r = (1, 0)
    for _ in range(n): r = gmul(r, x)
    return r
UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def gdivexact(x, y):
    """x/y if y | x in Z[i], else None."""
    n = gnorm(y)
    num = gmul(x, gconj(y))
    if num[0] % n or num[1] % n: return None
    return (num[0] // n, num[1] // n)

def gdivmod(x, y):
    n = gnorm(y)
    num = gmul(x, gconj(y))
    q = (round_div(num[0], n), round_div(num[1], n))
    r = gsub(x, gmul(q, y))
    return q, r

def round_div(a, b):
    # nearest integer to a/b (b>0)
    return (2 * a + b) // (2 * b)

def ggcd(x, y):
    while y != (0, 0):
        _, r = gdivmod(x, y)
        x, y = y, r
    return x

def gval(x, pi):
    """valuation of x at Gaussian prime pi (x != 0)."""
    v = 0
    while True:
        q = gdivexact(x, pi)
        if q is None: return v
        x = q; v += 1

def is_conj_primitive(x):
    """no rational integer >1 divides x and (1+i) does not divide x"""
    return math.gcd(x[0], x[1]) == 1 and (x[0] + x[1]) % 2 == 1

# ---------- primes ----------
def is_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def split_primes(lo, hi):
    return [p for p in range(lo, hi) if p % 4 == 1 and is_prime(p)]

def gaussian_prime_above(p):
    """pi = a+bi, a>b>0, a^2+b^2=p (p = 1 mod 4)."""
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2 and a > b > 0: return (a, b)
    raise ValueError(p)

# ---------- circles ----------
class Circle:
    """N = prod p_k^{e_k}, all p_k = 1 mod 4 (primitive circle: global gcd of all points is a unit)."""
    def __init__(self, primes, exps=None):
        self.primes = list(primes)
        self.exps = list(exps) if exps else [1] * len(self.primes)
        self.pis = [gaussian_prime_above(p) for p in self.primes]
        self.N = 1
        for p, e in zip(self.primes, self.exps): self.N *= p ** e
        self.W = math.log(self.N)
    def point(self, avec, unit=0):
        z = UNITS[unit]
        for pi, a, e in zip(self.pis, avec, self.exps):
            z = gmul(z, gmul(gpow(pi, a), gpow(gconj(pi), e - a)))
        return z
    def all_points(self):
        import itertools
        pts = []
        for avec in itertools.product(*[range(e + 1) for e in self.exps]):
            for u in range(4):
                pts.append((self.point(avec, u), avec, u))
        return pts
    def dist(self, a1, a2):
        return sum(abs(x - y) * math.log(p) for x, y, p in zip(a1, a2, self.primes))

def angle(z): return math.atan2(z[1], z[0]) % (2 * math.pi)

def best_clusters(circ, M):
    """for each window of M angularly consecutive points, the normalized chord constant
    C = |z_first - z_last| / N^(1/4); returns sorted list (C, window) best first."""
    pts = circ.all_points()
    pts.sort(key=lambda t: angle(t[0]))
    n = len(pts)
    out = []
    for i in range(n):
        win = [pts[(i + j) % n] for j in range(M)]
        ch2 = gnorm(gsub(win[0][0], win[-1][0]))
        out.append((math.sqrt(ch2) / circ.N ** 0.25, win))
    out.sort(key=lambda t: t[0])
    return out

def det_gauss(Mat):
    """exact determinant of a square matrix of Gaussian integers (fraction-free Bareiss over Z[i])."""
    n = len(Mat)
    A = [row[:] for row in Mat]
    sign = 1
    prev = (1, 0)
    for k in range(n - 1):
        if A[k][k] == (0, 0):
            for r in range(k + 1, n):
                if A[r][k] != (0, 0):
                    A[k], A[r] = A[r], A[k]; sign = -sign; break
            else:
                return (0, 0)
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = gsub(gmul(A[i][j], A[k][k]), gmul(A[i][k], A[k][j]))
                q = gdivexact(num, prev)
                assert q is not None
                A[i][j] = q
        prev = A[k][k]
    d = A[n - 1][n - 1]
    return d if sign == 1 else (-d[0], -d[1])
