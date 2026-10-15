# reflib.py -- referee's own exact helpers (written from the definitions in sharp.md Def 1.1 / L4K.md sec 1;
# nothing imported from the author's rigcompute_checks or the repository).
import math, itertools
from sympy import isprime, primerange, factorint

def m_of(q):  # |T_q| = q - chi_{-4}(q)
    return q - 1 if q % 4 == 1 else q + 1

def gprime(p):
    """pi = a+bi, a odd, b even, both > 0, a^2+b^2 = p."""
    for a in range(1, math.isqrt(p) + 1, 2):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2 and b > 0 and b % 2 == 0:
            return (a, b)
    raise ValueError(p)

def cmul(x, y, Q=None):
    a = x[0] * y[0] - x[1] * y[1]; b = x[0] * y[1] + x[1] * y[0]
    if Q: a %= Q; b %= Q
    return (a, b)

def Zv(pis, v):
    z = (1, 0)
    for (a, b), t in zip(pis, v):
        f = (a, b) if t > 0 else (a, -b)
        for _ in range(abs(t)):
            z = cmul(z, f)
    return z

def Dvals(z):  # s -> D_s ; s=0:Y (ratio 1), 1:X-Y (ratio i), 2:X (ratio -1), 3:X+Y (ratio -i)
    X, Y = z
    return [Y, X - Y, X, X + Y]

def val(x, q):
    if x == 0: return 10 ** 9
    x = abs(x); k = 0
    while x % q == 0: x //= q; k += 1
    return k

def norm_v(ps, v):
    N = 1
    for p, t in zip(ps, v): N *= p ** abs(t)
    return N

def capv(ps, v, s):
    return norm_v(ps, v) * (2 if s % 2 else 1)

def monos(r, K, half=False):
    out = []
    for v in itertools.product(range(-K, K + 1), repeat=r):
        k = sum(abs(t) for t in v)
        if 1 <= k <= K:
            if half:
                nz = [t for t in v if t]
                if nz[0] < 0: continue
            out.append(v)
    return out

def legendre(a, q):
    t = pow(a % q, (q - 1) // 2, q)
    return -1 if t == q - 1 else t

def nu_par(ps, q):
    return [1 if legendre(p, q) == -1 else 0 for p in ps]

class TorusLog:
    """Norm-one torus of (Z[i]/q)^*, as z -> z/conj z.  Represented concretely:
    split q: Z[i]/q = F_q via i -> r (r^2=-1); z/conj z -> (a+br)/(a-br) in F_q^*, cyclic order q-1.
    inert q: Z[i]/q = F_{q^2}; norm-one subgroup cyclic order q+1.  We pick a generator g of T and log by table.
    Returns logs relative to a generator; iota = log(i/conj i)?? NO: iota = log of i viewed in T, i.e. i itself
    (i has norm 1).  Collision of v with unit s at level 1: Z_v/conj Z_v = i^s."""
    def __init__(self, q):
        self.q = q; self.m = m = m_of(q)
        # brute-force elements of T as Gaussian residues: z/conj z for z = a+bi
        self.elems = {}
        fac = list(factorint(m).keys())
        # find generator by trying x = z/conj(z)
        def tor(a, b):
            nn = (a * a + b * b) % q
            inv = pow(nn, -1, q)
            return (((a * a - b * b) * inv) % q, (2 * a * b * inv) % q)
        def pw(x, e):
            r_ = (1, 0)
            while e:
                if e & 1: r_ = cmul(r_, x, q)
                x = cmul(x, x, q); e >>= 1
            return r_
        g = None
        for a in range(0, q):
            for b in range(1, q):
                if (a * a + b * b) % q == 0: continue
                t = tor(a, b)
                if all(pw(t, m // f) != (1, 0) for f in fac):
                    g = t; break
            if g: break
        self.g = g
        tab = {}; x = (1, 0)
        for k in range(m):
            tab[x] = k; x = cmul(x, g, q)
        assert len(tab) == m
        self.tab = tab; self.tor = tor
        self.iota = tab[(0, 1)]
        assert self.iota in (m // 4, 3 * m // 4)
    def L(self, z):  # log of z/conj z
        return self.tab[self.tor(z[0] % self.q, z[1] % self.q)]

def genuine_L(ps, q):
    T = TorusLog(q)
    return [T.L(gprime(p)) for p in ps], T.iota, T.m

def units(n):
    return [u for u in range(n) if math.gcd(u, n) == 1] if n > 1 else [0]

def eps_list(r):
    # sign vectors up to global sign: eps_0 = +1
    return [tuple([1] + [(-1 if (b >> j) & 1 else 1) for j in range(r - 1)]) for b in range(2 ** (r - 1))]

def class_mask(c, L, n, r):
    """bitmask over eps_list(r): bit k set iff c = u * eps_k * L mod n for some unit u (mod n)."""
    if n == 1: return (1 << (2 ** (r - 1))) - 1
    cs = tuple(x % n for x in c)
    mask = 0
    for k, e in enumerate(eps_list(r)):
        base = [(ei * x) % n for ei, x in zip(e, L)]
        for u in units(n):
            if tuple((u * x) % n for x in base) == cs:
                mask |= 1 << k; break
    return mask

from fractions import Fraction as Fr
def angle_ok_exact(ms, N):
    """Exact: exists x with |sin(x - s pi/4)| >= nu_s m_s / sqrt(N) for s=0..3 (nu_0=nu_2=1, nu_1=nu_3=2^-1/2).
    With y=2x: (cos y, sin y) in box [2m2^2/N-1, 1-2m0^2/N] x [m3^2/N-1, 1-m1^2/N]."""
    m0, m1, m2, m3 = ms
    xl, xh = Fr(2 * m2 * m2, N) - 1, 1 - Fr(2 * m0 * m0, N)
    yl, yh = Fr(m3 * m3, N) - 1, 1 - Fr(m1 * m1, N)
    if xl > xh or yl > yh: return False
    def cl(lo, hi): return 0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
    mn = cl(xl, xh) ** 2 + cl(yl, yh) ** 2
    mx = max(xl * xl, xh * xh) + max(yl * yl, yh * yh)
    return mn <= 1 <= mx
