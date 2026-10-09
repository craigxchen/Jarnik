"""Exact Gaussian-integer toolkit for round4/outside.md (pure Python integers).

A Gaussian integer is a tuple (a, b) = a + b i.
"""
from math import isqrt, atan2, pi, sqrt


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def conj(z):
    return (z[0], -z[1])


def sub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def _rdiv(a, n):
    return (2 * a + n) // (2 * n)


def divmod_g(z, w):
    n = norm(w)
    num = mul(z, conj(w))
    q = (_rdiv(num[0], n), _rdiv(num[1], n))
    r = sub(z, mul(q, w))
    assert norm(r) < n
    return q, r


def gcd_g(z, w):
    while w != (0, 0):
        _, r = divmod_g(z, w)
        z, w = w, r
    return z


def divides(w, z):
    n = norm(w)
    num = mul(z, conj(w))
    return num[0] % n == 0 and num[1] % n == 0


def exact_div(z, w):
    n = norm(w)
    num = mul(z, conj(w))
    assert num[0] % n == 0 and num[1] % n == 0
    return (num[0] // n, num[1] // n)


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def sqrt_m1(p):
    """x with x^2 = -1 mod p, p = 1 mod 4 prime."""
    for a in range(2, p):
        x = pow(a, (p - 1) // 4, p)
        if x * x % p == p - 1:
            return x
    raise ValueError


def gaussian_prime_above(p):
    """pi = a + b i with a^2 + b^2 = p, a > b > 0 (p = 1 mod 4)."""
    x = sqrt_m1(p)
    a, b = p, x
    while b * b > p:
        a, b = b, a % b
    c = isqrt(p - b * b)
    assert b * b + c * c == p
    u, v = max(b, c), min(b, c)
    return (u, v)


def gpow(z, e):
    r = (1, 0)
    for _ in range(e):
        r = mul(r, z)
    return r


def reps(N):
    """All (x, y) in Z^2 with x^2 + y^2 = N (exact, by factorisation)."""
    f = factor(N)
    base = [(1, 0)]
    for p, e in f.items():
        if p == 2:
            new = [mul(z, gpow((1, 1), e)) for z in base]
        elif p % 4 == 3:
            if e % 2:
                return []
            new = [(z[0] * p ** (e // 2), z[1] * p ** (e // 2)) for z in base]
        else:
            pi_ = gaussian_prime_above(p)
            new = []
            for z in base:
                for a in range(e + 1):
                    new.append(mul(z, mul(gpow(pi_, a), gpow(conj(pi_), e - a))))
        base = new
    out = set()
    for z in base:
        w = z
        for _ in range(4):
            out.add(w)
            w = mul(w, (0, 1))
    pts = sorted(out, key=lambda z: atan2(z[1], z[0]) % (2 * pi))
    for z in pts:
        assert norm(z) == N
    return pts


def clusters(N, C, minsize=3):
    """Maximal runs of consecutive points (by angle) whose angular span is <= C N^{-1/4}
    (arc length <= C sqrt(R)).  Floating point is used only to select runs."""
    pts = reps(N)
    m = len(pts)
    if m < minsize:
        return []
    ang = [atan2(z[1], z[0]) % (2 * pi) for z in pts]
    lim = C * N ** (-0.25)
    out = []
    j = 0
    for i in range(m):
        # extend window starting at i
        k = i
        while k + 1 < i + m:
            d = (ang[(k + 1) % m] - ang[i]) % (2 * pi)
            if d <= lim:
                k += 1
            else:
                break
        if k - i + 1 >= minsize:
            out.append([pts[t % m] for t in range(i, k + 1)])
    # keep maximal ones
    res = []
    for c in out:
        s = set(c)
        if not any(s < set(d) for d in out):
            if s not in [set(r) for r in res]:
                res.append(c)
    return res
