"""Exact arithmetic in T_q (norm-one units of Z[i]/q), q odd prime. Discrete logs w.r.t. a fixed generator."""
from sympy import primerange, isprime, sqrt_mod, primitive_root, legendre_symbol
import math

class Tq:
    def __init__(self, q):
        self.q = q
        self.split = (q % 4 == 1)
        self.m = q - 1 if self.split else q + 1
        q_ = q
        if self.split:
            self.c = int(sqrt_mod(q_ - 1, q_))  # c^2 = -1
            g = primitive_root(q_)
            self.log = {}
            x = 1
            for e in range(q_ - 1):
                self.log[x] = e; x = x * g % q_
        else:
            # find generator of norm-one group of F_{q^2}=F_q[i]
            m = self.m
            fac = factor_small(m)
            elems = []
            # z/zbar for z = a+bi ; enumerate a=1.., b=1
            for a in range(0, q_):
                for b in range(1, q_):
                    t = self.div_conj(a, b)
                    if t is None: continue
                    if all(self.pow(t, m // p) != (1, 0) for p in fac):
                        self.gen = t; break
                else: continue
                break
            self.log = {}
            x = (1, 0)
            for e in range(m):
                self.log[x] = e; x = self.mul(x, self.gen)
            assert len(self.log) == m
    def mul(self, x, y):
        q = self.q; return ((x[0]*y[0] - x[1]*y[1]) % q, (x[0]*y[1] + x[1]*y[0]) % q)
    def pow(self, x, e):
        r = (1, 0)
        while e:
            if e & 1: r = self.mul(r, x)
            x = self.mul(x, x); e >>= 1
        return r
    def div_conj(self, a, b):
        q = self.q; n = (a*a + b*b) % q
        if n == 0: return None
        inv = pow(n, -1, q)
        # (a+bi)/(a-bi) = (a+bi)^2/n
        return ((a*a - b*b) * inv % q, (2*a*b) * inv % q)
    def dlog_ratio(self, a, b):
        """log of (a+bi)/conj(a+bi) in Z/m, or None if not a unit."""
        q = self.q
        if self.split:
            x = (a + b*self.c) % q; y = (a - b*self.c) % q
            if x == 0 or y == 0: return None
            return self.log[x * pow(y, -1, q) % q]
        t = self.div_conj(a % q, b % q)
        if t is None: return None
        return self.log[t]
    def dlog_i(self):
        return self.dlog_ratio(1, 1)  # (1+i)/(1-i) = i

def factor_small(n):
    fs = []; p = 2
    while p * p <= n:
        if n % p == 0:
            fs.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: fs.append(n)
    return fs

def gauss_prime(p):
    """a+bi with a^2+b^2=p, p=1 mod 4 or p=2."""
    if p == 2: return (1, 1)
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)
    raise ValueError
