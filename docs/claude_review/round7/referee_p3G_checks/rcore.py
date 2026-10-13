"""Referee's independent core for p3G.md (no import of the author's code).

* exact Gaussian primes, tight relations of weight <= 2 (singles e_j, doubles 2e_j, pairs e_i +- e_j);
* tightness: |D_s(A_v)| divides L'_X (exact big-int modulus);
* orientation sigma by the exact integer formula of Lemma 2.3(2) (re-derived by the referee),
  sin^2 d = nu_s^2 D_s^2 / Norm(A) as an exact Fraction;
* closure with k <= 2 certificates via a quotient-group representation (signed components with
  status Z / Z2 / 0), written independently;
* an independent exact checker: sparse incremental echelon (HNF-type) lattice over Z.
"""
import math, random
from fractions import Fraction

MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

def is_prime(n):
    if n < 2: return False
    for q in MR_BASES:
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def gauss(p):
    """a + b i with a^2 + b^2 = p, a odd > 0, b even (sign of b as found)."""
    assert p % 4 == 1 and is_prime(p)
    c = 2
    while pow(c, (p - 1) // 2, p) != p - 1: c += 1
    x = pow(c, (p - 1) // 4, p)
    r0, r1 = p, x
    while r1 * r1 > p: r0, r1 = r1, r0 % r1
    a = r1; b = math.isqrt(p - a * a)
    assert a * a + b * b == p
    if a % 2 == 0: a, b = b, a
    return (a, b)

def odd_primes_upto(X):
    return [q for q in range(3, X + 1, 2) if all(q % t for t in range(3, math.isqrt(q) + 1, 2))]

def Lprime(X, ps):
    L = 1
    for q in odd_primes_upto(X):
        if q in ps: continue
        a = 1
        while q ** (a + 1) <= X: a += 1
        L *= q ** a
    return L

def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

def conj(z):
    return (z[0], -z[1])

def Ds(A):
    X, Y = A
    return (Y, X - Y, X, X + Y)

def sigma_exact(A, s):
    X, Y = A
    v = (X * Y, Y * Y - X * X, -X * Y, X * X - Y * Y)[s]
    return (v > 0) - (v < 0)

def sin2d(A, s):
    X, Y = A
    D = Ds(A)[s]
    return Fraction(D * D, X * X + Y * Y) * (Fraction(1, 2) if s % 2 else 1)

def tight_list(gp, L):
    """All tight (v, s) of weight <= 2. v as tuple of (col, coef). Returns list of dicts."""
    r = len(gp); out = []
    def test(v, A):
        for s, D in enumerate(Ds(A)):
            if D != 0 and L % abs(D) == 0:
                out.append((v, s, sigma_exact(A, s), sin2d(A, s)))
    for j in range(r):
        test(((j, 1),), gp[j])
        test(((j, 2),), gmul(gp[j], gp[j]))
    for i in range(r):
        zi = gp[i]
        for j in range(i + 1, r):
            test(((i, 1), (j, 1)), gmul(zi, gp[j]))
            test(((i, 1), (j, -1)), gmul(zi, conj(gp[j])))
    return out

def rand_primes_like_author(r, lo, hi, rng):
    """Same RNG consumption as the description in p3G.md 4.2 (randrange | 1, p = 1 mod 4)."""
    out = set()
    while len(out) < r:
        n = rng.randrange(lo, hi) | 1
        if n % 4 == 1 and is_prime(n): out.add(n)
    return sorted(out)
