"""Minimal exact Gaussian-integer helpers (pure integers, no sympy).

A Gaussian integer is a tuple (a, b) meaning a + b i.
"""
from math import isqrt


def add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def sub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def _round_div(a, n):
    # nearest integer to a/n, n > 0
    return (2 * a + n) // (2 * n)


def divmod_g(z, w):
    """Euclidean division z = q w + r with N(r) < N(w)."""
    n = norm(w)
    if n == 0:
        raise ZeroDivisionError
    num = mul(z, conj(w))
    q = (_round_div(num[0], n), _round_div(num[1], n))
    r = sub(z, mul(q, w))
    assert norm(r) < n
    return q, r


def divides(w, z):
    """True iff w | z in Z[i] (w != 0)."""
    n = norm(w)
    num = mul(z, conj(w))
    return num[0] % n == 0 and num[1] % n == 0


def exact_div(z, w):
    n = norm(w)
    num = mul(z, conj(w))
    assert num[0] % n == 0 and num[1] % n == 0, (z, w)
    return (num[0] // n, num[1] // n)


def gcd_g(z, w):
    while w != (0, 0):
        _, r = divmod_g(z, w)
        z, w = w, r
    return z


def is_unit(z):
    return norm(z) == 1


def is_prime_int(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return False
        f += 2
    return True


def factor_int(n):
    """Trial-division factorisation of a positive integer: dict p -> e."""
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def split_prime_factor(p):
    """For a prime p = 1 mod 4 return a Gaussian prime pi with N(pi) = p."""
    assert p % 4 == 1
    # find t with t^2 = -1 mod p
    for a in range(2, p):
        t = pow(a, (p - 1) // 4, p)
        if (t * t) % p == p - 1:
            break
    g = gcd_g((p, 0), (t, 1))
    assert norm(g) == p
    return g


def valuation(pi, z):
    """Largest a with pi^a | z (z != 0)."""
    a = 0
    while divides(pi, z):
        z = exact_div(z, pi)
        a += 1
    return a


def gpow(z, k):
    out = (1, 0)
    for _ in range(k):
        out = mul(out, z)
    return out


def lattice_points(N):
    """All Gaussian integers of norm N (brute force)."""
    pts = []
    r = isqrt(N)
    for x in range(-r, r + 1):
        y2 = N - x * x
        y = isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    return pts
