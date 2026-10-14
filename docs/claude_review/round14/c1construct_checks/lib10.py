# Shared exact helpers for the L4(K) tests.
# Labels at prime q (not dividing N): ell in Z/m, m = q - chi(q), chi = chi_{-4}; n = m/4.
# v collides (some unit) iff <v,ell> = 0 mod n; unit s = (<v,ell> mod m)/n.
# Parity-respecting: ell_j = nu_j(q) mod 2, nu_j(q) = [(p_j/q) = -1].
from sympy import primerange, isprime, legendre_symbol
import itertools, math

def mq(q):
    return q - (1 if q % 4 == 1 else -1)

def nu(ps, q):
    return [1 if legendre_symbol(p % q, q) == -1 else 0 for p in ps]

def split_primes(lo, hi):
    return [p for p in primerange(lo, hi + 1) if p % 4 == 1]

def monos(r, K, half=True):
    """all v in Z^r, 1<=||v||_1<=K; if half, one of each +-pair (first nonzero coord > 0)."""
    out = []
    def rec(i, rem, cur):
        if i == r:
            if any(cur):
                if not half or next(t for t in cur if t) > 0:
                    out.append(tuple(cur))
            return
        for t in range(-rem, rem + 1):
            cur.append(t); rec(i + 1, rem - abs(t), cur); cur.pop()
    rec(0, K, [])
    return out

def gauss_prime(p):
    """a+bi with a^2+b^2=p, a odd >0, b even >0 (fixed normalisation)."""
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2 and a % 2 == 1:
            return (a, b)
    raise ValueError(p)

def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gconj(x):
    return (x[0], -x[1])

def Zv(pis, v):
    z = (1, 0)
    for pi, t in zip(pis, v):
        f = pi if t > 0 else gconj(pi)
        for _ in range(abs(t)):
            z = gmul(z, f)
    return z

def Ds(z, s):
    X, Y = z
    return [Y, X - Y, X, X + Y][s % 4]

def wv(ps, v):
    return sum(abs(t) * math.log(p) for p, t in zip(ps, v))

def budget(ps, v, s):
    return wv(ps, v) / 2 + (0.5 * math.log(2) if s % 2 else 0.0)

# --- discrete log on T_q (norm-one torus of (Z[i]/q)^*) ---
class Torus:
    """T_q for odd prime q: elements z/conj z; cyclic of order m. log: T_q -> Z/m."""
    def __init__(self, q):
        self.q = q; self.m = mq(q)
        # find generator: brute force over Gaussian residues
        self.table = {}
        q = self.q; m = self.m
        fac = [f for f in range(2, m + 1) if m % f == 0 and isprime(f)]
        for a in range(q):
            for b in range(q):
                n = (a * a + b * b) % q
                if n == 0: continue
                inv = pow(n, -1, q)
                # t = z / conj z = z^2 / N(z)
                t = (((a * a - b * b) * inv) % q, (2 * a * b * inv) % q)
                if all(self.pw(t, m // f) != (1, 0) for f in fac):
                    self.g = t
                    break
            else:
                continue
            break
        x = (1, 0)
        for k in range(m):
            self.table[x] = k
            x = self.mul(x, self.g)
        assert len(self.table) == m
        self.iota = self.table[(0, 1)]
        assert self.iota % (m // 4) == 0
    def mul(self, x, y):
        q = self.q
        return ((x[0] * y[0] - x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)
    def pw(self, x, e):
        r = (1, 0)
        while e:
            if e & 1: r = self.mul(r, x)
            x = self.mul(x, x); e >>= 1
        return r
    def L(self, z):
        """log(z / conj z) in Z/m, z Gaussian integer prime to q (as tuple)."""
        q = self.q; a, b = z[0] % q, z[1] % q
        n = (a * a + b * b) % q
        assert n != 0
        inv = pow(n, -1, q)
        t = (((a * a - b * b) * inv) % q, (2 * a * b * inv) % q)
        return self.table[t]
