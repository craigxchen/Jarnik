"""Small exact Gaussian-integer toolkit for the referee-2 checks (Python ints only)."""
from math import isqrt


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
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


def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def power(z, k):
    r = (1, 0)
    for _ in range(k):
        r = mul(r, z)
    return r


def divides(a, z):
    n = norm(a); w = mul(z, conj(a))
    return w[0] % n == 0 and w[1] % n == 0


def div(z, a):
    n = norm(a); w = mul(z, conj(a))
    assert w[0] % n == 0 and w[1] % n == 0
    return (w[0] // n, w[1] // n)


def gaussian_prime_over(p):
    """primary Gaussian prime a+bi over p = 1 mod 4: a odd, b even, a+b = 1 mod 4"""
    for x in range(1, isqrt(p) + 1):
        y2 = p - x * x
        y = isqrt(y2)
        if y * y == y2:
            for (a, b) in ((x, y), (y, x), (-x, y), (x, -y), (-x, -y), (y, -x), (-y, x), (-y, -x)):
                if a % 2 == 1 and b % 2 == 0 and (a + b) % 4 == 1:
                    return (a, b)
    raise ValueError


def mulmod(z, w, m):
    return ((z[0] * w[0] - z[1] * w[1]) % m, (z[0] * w[1] + z[1] * w[0]) % m)


def invmod(z, m):
    n = norm(z) % m
    ni = pow(n, -1, m)
    c = conj(z)
    return ((c[0] * ni) % m, (c[1] * ni) % m)


def powmod(z, k, m):
    if k < 0:
        z = invmod(z, m); k = -k
    r = (1, 0); b = (z[0] % m, z[1] % m)
    while k:
        if k & 1:
            r = mulmod(r, b, m)
        b = mulmod(b, b, m); k >>= 1
    return r


def lattice_points(N):
    pts = []
    for x in range(-isqrt(N), isqrt(N) + 1):
        y2 = N - x * x
        y = isqrt(y2)
        if y * y == y2:
            pts.append((x, y))
            if y:
                pts.append((x, -y))
    return pts


def Ksize(q, a):
    """order of the norm-one subgroup K of (Z[i]/q^a)^*, q odd"""
    return (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)


def Kclass(t, q, a):
    """class of a norm-one element t in K/K^2 = {+1,-1}: t^(|K|/2) mod q^a is (1,0) or (-1,0)"""
    m = q ** a
    v = powmod(t, Ksize(q, a) // 2, m)
    if v == (1, 0):
        return 1
    if v == ((-1) % m, 0):
        return -1
    raise AssertionError("not in K", t, q, a, v)
