"""Referee toolkit (written independently of stability_checks/stab_lib.py).

Gaussian integers are Python complex-free pairs (a, b) = a + b i with exact integer arithmetic.
F*(s) of two_thirds Theorem 5.2 is returned as an exact integer exponent vector
    8 F*(s) = sum_q n_q log q        (n_q integers),
so every inequality "... + 8F*(s) <= ..." can be exponentiated and checked in integers.
The split-prime minimum Psi_p is computed by a shortest-path DP over the cut values S_1 < ... < S_e
(different code path from stab_lib.psi8), and cross-checked against brute force over compositions.
"""
from math import comb, isqrt
from itertools import combinations


# ------------------------------------------------------------------ Gaussian integers
def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gsub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def gconj(z):
    return (z[0], -z[1])


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def gdivexact(z, w):
    n = gnorm(w)
    t = gmul(z, gconj(w))
    if t[0] % n or t[1] % n:
        return None
    return (t[0] // n, t[1] // n)


def gmod(z, w):
    n = gnorm(w)
    t = gmul(z, gconj(w))
    # nearest-integer quotient
    q = ((2 * t[0] + n) // (2 * n), (2 * t[1] + n) // (2 * n))
    return gsub(z, gmul(q, w))


def ggcd(z, w):
    while w != (0, 0):
        z, w = w, gmod(z, w)
    return z


def ggcd_list(zs):
    g = zs[0]
    for z in zs[1:]:
        g = ggcd(g, z)
        if gnorm(g) == 1:
            return g
    return g


def gpow(z, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, z)
    return r


# ------------------------------------------------------------------ rational integers
def factorize(n):
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


def sieve(n):
    b = bytearray([1]) * (n + 1)
    b[0] = b[1] = 0
    for i in range(2, isqrt(n) + 1):
        if b[i]:
            b[i * i::i] = bytearray(len(b[i * i::i]))
    return [i for i in range(n + 1) if b[i]]


_PI = {}


def split_prime(p):
    """A Gaussian prime of norm p (p = 1 mod 4), via Cornacchia-type descent."""
    if p in _PI:
        return _PI[p]
    for c in range(2, p):
        if pow(c, (p - 1) // 2, p) == p - 1:
            t = pow(c, (p - 1) // 4, p)
            break
    a, b = p, t
    while b * b > p:
        a, b = b, a % b
    x = b
    y = isqrt(p - x * x)
    assert x * x + y * y == p
    _PI[p] = (x, y)
    return (x, y)


def lattice_points(fac):
    """All Gaussian integers with norm prod q^e (fac), as a list (unsorted)."""
    base = [(1, 0)]
    for q, e in fac.items():
        if q == 2:
            base = [gmul(z, gpow((1, 1), e)) for z in base]
        elif q % 4 == 3:
            if e % 2:
                return []
            base = [gmul(z, (q ** (e // 2), 0)) for z in base]
        else:
            pi = split_prime(q)
            pb = gconj(pi)
            base = [gmul(z, gmul(gpow(pi, a), gpow(pb, e - a))) for z in base for a in range(e + 1)]
    out = set()
    for z in base:
        w = z
        for _ in range(4):
            out.add(w)
            w = (-w[1], w[0])
    return list(out)


def vpi(z, pi):
    """valuation of Gaussian integer z at Gaussian prime pi."""
    a = 0
    while True:
        w = gdivexact(z, pi)
        if w is None:
            return a
        z = w
        a += 1


# ------------------------------------------------------------------ F*(s), exact exponent vector
def E(n, Q):
    if n <= 1:
        return 0
    b, r = divmod(n, Q)
    return Q * comb(b, 2) + r * b


def depth(n, base, q):
    """sum_{a>=1} E(n, base*q^(a-1))."""
    t, Q = 0, base
    while Q < n:
        t += E(n, Q)
        Q *= q
    return t


def psi8_dp(M, p):
    """8 * Psi_p(M): min over e>=1 and level sizes n_0..n_e >= 1 (sum M) of
       8*sum_t depth(n_t) + sum_tau (2 S_tau - M)^2.  DP over the last cut value S."""
    w = [8 * depth(n, p - 1, p) for n in range(M + 1)]
    INF = float('inf')
    best = [INF] * (M + 1)          # best[S]: cost of levels below cut S, cuts counted, last cut at S
    for S in range(1, M):
        best[S] = w[S] + (2 * S - M) ** 2
    for S in range(1, M):
        if best[S] == INF:
            continue
        for T in range(S + 1, M):
            c = best[S] + w[T - S] + (2 * T - M) ** 2
            if c < best[T]:
                best[T] = c
    return min(best[S] + w[M - S] for S in range(1, M))


def psi8_brute(M, p):
    """Same, brute force over all compositions of M into >= 2 positive parts."""
    w = [8 * depth(n, p - 1, p) for n in range(M + 1)]
    best = None
    for k in range(1, M):
        for cuts in combinations(range(1, M), k):
            parts = [cuts[0]] + [cuts[i] - cuts[i - 1] for i in range(1, k)] + [M - cuts[-1]]
            v = sum(w[x] for x in parts) + sum((2 * S - M) ** 2 for S in cuts)
            if best is None or v < best:
                best = v
    return best


_F8 = {}


def F8vec(s, brute=False):
    """dict q -> n_q with 8F*(s) = sum n_q log q (exact integers)."""
    key = (s, brute)
    if key in _F8:
        return _F8[key]
    vec = {2: 4 * comb(s, 2)}
    for q in sieve(s + 2):
        if q == 2:
            continue
        if q % 4 == 3:
            v = 8 * depth(s, q + 1, q)
        else:
            free = 8 * depth(s, q - 1, q)
            if q - 1 >= s:
                v = 0
            else:
                ps = psi8_brute(s, q) if brute else psi8_dp(s, q)
                v = min(free, ps)
        if v:
            vec[q] = v
    _F8[key] = vec
    return vec


def F8int(s):
    """exp(8 F*(s)) as an integer."""
    r = 1
    for q, n in F8vec(s).items():
        r *= q ** n
    return r


def F8float(s):
    from math import log
    return sum(n * log(q) for q, n in F8vec(s).items())
