"""Shared helpers for round6/sharp.md: Gaussian integers modulo a prime power Q = q^a (q odd),
the norm-one group T_Q, discrete logs, square roots.  Pure Python, exact integer arithmetic.

A residue is a pair (x, y) meaning x + i y in Z[i]/Q."""


def is_prime(n):
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


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def chi4(q):
    return 1 if q % 4 == 1 else -1


def m_of(q, a):
    """order of the norm-one group T_{q^a} = (q - chi(q)) q^(a-1)"""
    return (q - chi4(q)) * q ** (a - 1)


def mul(u, v, Q):
    return ((u[0] * v[0] - u[1] * v[1]) % Q, (u[0] * v[1] + u[1] * v[0]) % Q)


def conj(u, Q):
    return (u[0] % Q, (-u[1]) % Q)


def norm(u, Q):
    return (u[0] * u[0] + u[1] * u[1]) % Q


def gpow(u, e, Q):
    r = (1, 0)
    b = (u[0] % Q, u[1] % Q)
    while e > 0:
        if e & 1:
            r = mul(r, b, Q)
        b = mul(b, b, Q)
        e >>= 1
    return r


def inv(u, Q):
    n = norm(u, Q)
    ni = pow(n, -1, Q)
    c = conj(u, Q)
    return ((c[0] * ni) % Q, (c[1] * ni) % Q)


def is_unit(u, q, Q):
    return norm(u, Q) % q != 0


def order_in_T(t, m, Q):
    """order of t in the cyclic group T_Q of order m"""
    o = m
    for p in factor(m):
        while o % p == 0 and gpow(t, o // p, Q) == (1, 0):
            o //= p
    return o


def generator_T(q, a, seed=0):
    """a generator of T_{q^a}: take t = r/conj(r) for units r until the order is m"""
    import random
    Q = q ** a
    m = m_of(q, a)
    rng = random.Random(seed)
    while True:
        r = (rng.randrange(Q), rng.randrange(Q))
        if not is_unit(r, q, Q):
            continue
        t = mul(r, inv(conj(r, Q), Q), Q)
        if order_in_T(t, m, Q) == m:
            return t


def dlog_table(g, m, Q):
    """dict t -> log_g t for all t in T_Q (only for small m)"""
    tab = {}
    t = (1, 0)
    for k in range(m):
        tab[t] = k
        t = mul(t, g, Q)
    assert t == (1, 0)
    return tab


def hilbert90(t, q, Q, rng):
    """r with r/conj(r) = t (t in T_Q): r = u + t*conj(u) for a u making r a unit"""
    while True:
        u = (rng.randrange(Q), rng.randrange(Q))
        r = ((u[0] + mul(t, conj(u, Q), Q)[0]) % Q, (u[1] + mul(t, conj(u, Q), Q)[1]) % Q)
        if is_unit(r, q, Q):
            assert mul(r, inv(conj(r, Q), Q), Q) == t
            return r


def sqrt_mod_prime_power(n, q, a):
    """a square root of the unit n modulo q^a (q odd), or None"""
    n %= q
    if n == 0:
        return None
    if pow(n, (q - 1) // 2, q) != 1:
        return None
    # Tonelli-Shanks is overkill for our sizes: brute force mod q, then Hensel
    s = next(x for x in range(1, q) if (x * x - n) % q == 0)
    return s


def sqrt_unit(n, q, a):
    Q = q ** a
    s0 = sqrt_mod_prime_power(n % q, q, 1)
    if s0 is None:
        return None
    s = s0
    mod = q
    for _ in range(1, a):
        mod *= q
        # Newton step: s <- s - (s^2 - n)/(2s)
        s = (s - (s * s - n) * pow(2 * s, -1, mod)) % mod
    assert (s * s - n) % Q == 0
    return s


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1
