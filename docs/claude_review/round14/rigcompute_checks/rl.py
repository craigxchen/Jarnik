# rl.py -- shared exact helpers for round 14 (rigcompute).  Independent re-implementation (no imports from the
# repository).  Conventions follow round10/L4K.md section 1.
import math, itertools, random, subprocess, os
import numpy as np
from sympy import primerange, isprime, legendre_symbol, factorint

HERE = os.path.dirname(os.path.abspath(__file__))

def mq(q):  # m_q = q - chi_{-4}(q)
    return q - (1 if q % 4 == 1 else -1)

def nu(ps, q):
    return [1 if legendre_symbol(p % q, q) == -1 else 0 for p in ps]

def split_primes(lo, hi):
    return [p for p in primerange(lo, hi + 1) if p % 4 == 1]

def monos_norm(r, K):
    """all v in Z^r, 1<=||v||_1<=K, one of each +-pair: LAST nonzero coordinate > 0."""
    out = []
    def rec(i, rem, cur):
        if i == r:
            if any(cur):
                last = [t for t in cur if t][-1]
                if last > 0: out.append(tuple(cur))
            return
        for t in range(-rem, rem + 1):
            cur.append(t); rec(i + 1, rem - abs(t), cur); cur.pop()
    rec(0, K, [])
    return out

def gauss_prime(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a * a; b = math.isqrt(b2)
        if b * b == b2 and a % 2 == 1:
            return (a, b)
    raise ValueError(p)

def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def Zv(pis, v):
    z = (1, 0)
    for pi, t in zip(pis, v):
        f = pi if t > 0 else (pi[0], -pi[1])
        for _ in range(abs(t)): z = gmul(z, f)
    return z

def Ds(z, s):
    X, Y = z
    return (Y, X - Y, X, X + Y)[s % 4]

def vq(x, q):
    x = abs(x); k = 0
    while x and x % q == 0: x //= q; k += 1
    return k

def cap(ps, v, s):
    """exact budget: collision product P allowed iff P^2 <= cap = Norm(Z_v) * (2 if s odd)."""
    N = 1
    for p, t in zip(ps, v): N *= p ** abs(t)
    return N * (2 if s % 2 else 1)

class Torus:
    """norm-one torus T_q = {z/conj z} in (Z[i]/q)^*, cyclic of order m.  dlog via table of a generator."""
    def __init__(self, q):
        self.q = q; m = self.m = mq(q)
        fac = list(factorint(m).keys())
        g = None
        for a in range(1, q):
            for b in range(1, q):
                nn = (a * a + b * b) % q
                if nn == 0: continue
                inv = pow(nn, -1, q)
                t = (((a * a - b * b) * inv) % q, (2 * a * b * inv) % q)
                if all(self.pw(t, m // f) != (1, 0) for f in fac):
                    g = t; break
            if g: break
        self.g = g; self.table = {}
        x = (1, 0)
        for k in range(m):
            self.table[x] = k; x = self.mul(x, g)
        assert len(self.table) == m
        self.iota = self.table[(0, 1)]
        assert self.iota in (m // 4, 3 * m // 4)
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
        q = self.q; a, b = z[0] % q, z[1] % q
        nn = (a * a + b * b) % q; assert nn
        inv = pow(nn, -1, q)
        return self.table[(((a * a - b * b) * inv) % q, (2 * a * b * inv) % q)]

def units_mod(n):
    return [u for u in range(1, n + 1) if math.gcd(u, n) == 1] if n > 1 else [1]

def genuine_class(Lc, n):
    """{u*L mod n : u unit mod n} as a set of tuples (level-1 genuine class)."""
    return set(tuple((u * x) % n for x in Lc) for u in units_mod(n))

def orbit_reps_U1(m, parity):
    """orbit representatives (with orbit sizes) of x in Z/m, x = parity mod 2, under U1 = {u unit mod m, u=1 mod 4}."""
    U1 = [u for u in range(1, m, 4) if math.gcd(u, m) == 1]
    seen = bytearray(m); reps = []
    for x in range(parity, m, 2):
        if seen[x]: continue
        orb = set((u * x) % m for u in U1)
        for y in orb: seen[y] = 1
        reps.append((x, len(orb)))
    return reps

def _mk(t, o, r, L):
    nmk = int(t[o]); mks = []; p = o + 1
    for _ in range(nmk):
        mks.append((int(t[p]), tuple(int(x) for x in t[p + 1:p + 1 + r]))); p += 1 + r
    return dict(mks=mks, el=tuple(int(x) for x in t[p:p + L]))

def run_csp(M, n, r, mode, maxwit, maxnodes, par, col0, V, forb, cls, timeout=None, keep_all=False, rec=None):
    """call ./csp ; V list of sign-normalised tuples, forb list of residue lists;
    cls: list of tuples mod n (mask 1) or dict tuple->mask.  mode 3 returns signatures (keep_all forced)."""
    if mode == 3: keep_all = True
    if not isinstance(cls, dict): cls = {c: 1 for c in cls}
    lines = [f"{M} {n} {r} {mode} {maxwit} {maxnodes}", " ".join(map(str, par)), str(len(col0))]
    lines += [f"{x} {w}" for x, w in col0]
    if rec is None: rec = [15] * len(V)
    keep = [(v, f, rc) for v, f, rc in zip(V, forb, rec) if f or keep_all]
    lines.append(str(len(keep)))
    lines += [" ".join(map(str, v)) + f" {len(f)} " + " ".join(map(str, f)) + (f" {rc}" if mode == 3 else "") for v, f, rc in keep]
    lines.append(str(len(cls)))
    lines += [" ".join(map(str, c)) + f" {mk}" for c, mk in sorted(cls.items())]
    out = subprocess.run([os.path.join(HERE, "csp")], input="\n".join(lines) + "\n", capture_output=True,
                         text=True, timeout=timeout).stdout
    wits = []; res = None; sigs = []
    for ln in out.splitlines():
        if ln.startswith("W "): wits.append(tuple(int(x) for x in ln.split()[1:]))
        elif ln.startswith("S "):
            t = ln.split(); L = int(t[5])
            sigs.append(dict(orm=int(t[1]), none=int(t[2]), notL=int(t[3]), cnt=int(t[4]),
                             rep=tuple(int(x) for x in t[6:6 + r]), repN=tuple(int(x) for x in t[6 + r:6 + 2 * r]),
                             **_mk(t, 6 + 2 * r, r, L)))
        elif ln.startswith("RESULT"):
            t = ln.split(); res = dict(status=t[1], nodes=int(t[2]), leaves=int(t[3]), nongen=int(t[4]), raw=int(t[5]))
    if mode == 3: return res, sigs
    return res, wits
