"""Referee's own exact arithmetic in T_q = norm-one units of Z[i]/q (q odd prime). Independent of the author's tgroup.py.
Labels: L(z) = discrete log of z/conj(z) w.r.t. a generator g of T_q; IOTA = L(1+i) = log(i)."""
import math
from sympy import primitive_root, sqrt_mod, factorint

class T:
    def __init__(self, q):
        self.q = q; self.m = q - 1 if q % 4 == 1 else q + 1
        if q % 4 == 1:
            self.c = int(sqrt_mod(q - 1, q)); g = primitive_root(q)
            self.tab = {}; x = 1
            for e in range(q - 1): self.tab[x] = e; x = x * g % q
        else:
            primes = list(factorint(self.m))
            gen = None
            for a in range(q):
                for b in range(1, q):
                    t = self._ratio(a, b)
                    if t is None: continue
                    if all(self._pw(t, self.m // p) != (1, 0) for p in primes): gen = t; break
                if gen: break
            self.tab = {}; x = (1, 0)
            for e in range(self.m): self.tab[x] = e; x = self._mul(x, gen)
        self.iota = self.L(1, 1)
        assert self.iota in (self.m // 4, 3 * self.m // 4)
    def _mul(self, x, y):
        q = self.q; return ((x[0]*y[0] - x[1]*y[1]) % q, (x[0]*y[1] + x[1]*y[0]) % q)
    def _pw(self, x, e):
        r = (1, 0)
        while e:
            if e & 1: r = self._mul(r, x)
            x = self._mul(x, x); e >>= 1
        return r
    def _ratio(self, a, b):  # (a+bi)/(a-bi) in F_q[i], q = 3 mod 4
        q = self.q; n = (a*a + b*b) % q
        if n == 0: return None
        inv = pow(n, -1, q)
        return ((a*a - b*b) * inv % q, 2*a*b*inv % q)
    def L(self, a, b):
        q = self.q
        if q % 4 == 1:
            x = (a + b*self.c) % q; y = (a - b*self.c) % q
            if x == 0 or y == 0: return None
            return self.tab[x * pow(y, -1, q) % q]
        t = self._ratio(a % q, b % q)
        return None if t is None else self.tab[t]
    def unit(self, val):
        """s in Z/4 with val = s*iota mod m, or None."""
        val %= self.m
        if val % (self.m // 4): return None
        u = val // (self.m // 4)            # val = u*m/4
        return u if self.iota == self.m // 4 else (3 * u) % 4

def two_sq(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return a, b

def D(a, b, s):
    return [b, a - b, a, a + b][s % 4]
