"""Exact Gaussian-integer helpers (Python ints) for round 5."""
from math import isqrt


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def two_squares(p):
    """p = 1 mod 4 prime -> (a, b) with a^2 + b^2 = p, a odd > 0, b even > 0 (Cornacchia)."""
    assert p % 4 == 1
    # square root of -1 mod p
    for g in range(2, p):
        t = pow(g, (p - 1) // 4, p)
        if t * t % p == p - 1:
            break
    a, b = p, t
    lim = isqrt(p)
    while b > lim:
        a, b = b, a % b
    x = b
    y = isqrt(p - x * x)
    assert x * x + y * y == p
    if x % 2 == 0:
        x, y = y, x
    return (x, y)


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gpow(z, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, z)
    return r


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def gmod(z, m):
    return (z[0] % m, z[1] % m)


def ginv_mod(z, m):
    """inverse of z in (Z[i]/m)^*, m odd: z^{-1} = conj(z) * N(z)^{-1}"""
    n = gnorm(z) % m
    ni = pow(n, -1, m)
    c = gconj(z)
    return ((c[0] * ni) % m, (c[1] * ni) % m)


def gmulmod(z, w, m):
    return ((z[0] * w[0] - z[1] * w[1]) % m, (z[0] * w[1] + z[1] * w[0]) % m)


def gpowmod(z, k, m):
    r = (1, 0)
    if k < 0:
        z = ginv_mod(z, m); k = -k
    b = gmod(z, m)
    while k:
        if k & 1:
            r = gmulmod(r, b, m)
        b = gmulmod(b, b, m)
        k >>= 1
    return r


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def odd_prime_powers_upto(Q0):
    out = []
    for q in primes_upto(Q0):
        if q == 2:
            continue
        a, qa = 1, q
        while qa <= Q0:
            out.append((q, a, qa))
            a += 1; qa *= q
    return out
