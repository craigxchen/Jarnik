"""Referee's own library (independent of rigtheory_checks). Genuine labels via explicit residue fields."""
import math
from sympy import isprime, factorint, primerange
def two_sq(p):
    for a in range(1, math.isqrt(p)+1):
        b2 = p-a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)
def sqrtm1(q):
    for x in range(2, q):
        if pow(x, (q-1)//2, q) == q-1: return pow(x, (q-1)//4, q)
class Tor:
    """norm-one torus mod q (prime level), cyclic of order m; elements represented as ints (split) or pairs (inert)."""
    def __init__(self, q):
        self.q = q; self.split = (q % 4 == 1)
        self.m = q-1 if self.split else q+1
        self.n = self.m//4
        fac = list(factorint(self.m))
        if self.split:
            self.I = sqrtm1(q)  # i -> I in first factor
            for g in range(2, q):
                if all(pow(g, self.m//f, q) != 1 for f in fac): break
            self.g = g
            # dlog table
            self.tab = {}; x = 1
            for k in range(self.m): self.tab[x] = k; x = x*g % q
            li = self.tab[self.I]
        else:
            # F_q[i]; torus elements (a,b) with a^2+b^2=1
            def mul(z, w): return ((z[0]*w[0]-z[1]*w[1]) % q, (z[0]*w[1]+z[1]*w[0]) % q)
            self.mul = mul
            pts = [(a, b) for a in range(q) for b in range(q) if (a*a+b*b) % q == 1]
            assert len(pts) == self.m
            for g in pts:
                x = (1, 0); ordr = 0; ok = True
                for f in fac:
                    y = (1, 0); e = self.m//f; base = g
                    # pow
                    r_ = (1, 0); b_ = g; ee = e
                    while ee:
                        if ee & 1: r_ = mul(r_, b_)
                        b_ = mul(b_, b_); ee >>= 1
                    if r_ == (1, 0): ok = False; break
                if ok: break
            self.g = g; self.tab = {}; x = (1, 0)
            for k in range(self.m): self.tab[x] = k; x = mul(x, g)
            li = self.tab[(0, 1)]
        assert li in (self.n, 3*self.n)
        self.sign = 1 if li == self.n else -1
    def ratio_log(self, z):
        """dlog of z/conj(z), normalised so log i = m/4"""
        q = self.q; a, b = z
        if self.split:
            num = (a + b*self.I) % q; den = (a - b*self.I) % q
            t = num*pow(den, -1, q) % q
            k = self.tab[t]
        else:
            nz = (a*a+b*b) % q; inv = pow(nz, -1, q)
            # z/conj z = z^2 / Norm z
            z2 = ((a*a-b*b) % q, (2*a*b) % q)
            t = (z2[0]*inv % q, z2[1]*inv % q)
            k = self.tab[t]
        return (self.sign*k) % self.m
def genuine(ps, q):
    T = Tor(q); return T.m, [T.ratio_log(two_sq(p)) for p in ps]
def nu(ps, q): return [1 if pow(p % q, (q-1)//2, q) == q-1 else 0 for p in ps]
def Zv(pis, v):
    X, Y = 1, 0
    for (a, b), e in zip(pis, v):
        if e < 0: b = -b
        for _ in range(abs(e)): X, Y = X*a-Y*b, X*b+Y*a
    return X, Y
def Ds(z, s): X, Y = z; return (Y, X-Y, X, X+Y)[s % 4]
def VK(r, K):
    out = []
    def rec(i, rem, cur):
        if i == r:
            if any(cur): out.append(tuple(cur))
            return
        for x in range(-rem, rem+1): rec(i+1, rem-abs(x), cur+[x])
    rec(0, K, []); return out
def vq(x, q):
    x = abs(x); e = 0
    while x % q == 0: x //= q; e += 1
    return e
