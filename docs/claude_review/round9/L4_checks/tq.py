"""Exact arithmetic in T_Q (norm-one units of Z[i]/Q), Q = q^a, q odd prime. Discrete logs by table."""
import math, random

def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def primes_upto(n):
    return [p for p in range(2, n+1) if is_prime(p)]

def chi4(q): return 1 if q % 4 == 1 else -1

def gmul(a, b, Q):
    return ((a[0]*b[0]-a[1]*b[1]) % Q, (a[0]*b[1]+a[1]*b[0]) % Q)

def gpow(a, e, Q):
    r = (1, 0)
    while e:
        if e & 1: r = gmul(r, a, Q)
        a = gmul(a, a, Q); e >>= 1
    return r

def factor(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d,0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n,0)+1
    return f

def ratio(r, Q):
    """r/conj(r) = r^2 / Norm(r) mod Q."""
    n = (r[0]*r[0]+r[1]*r[1]) % Q
    inv = pow(n, -1, Q)
    r2 = gmul(r, r, Q)
    return (r2[0]*inv % Q, r2[1]*inv % Q)

class TQ:
    def __init__(self, q, a):
        self.q, self.a, self.Q = q, a, q**a
        self.m = (q - chi4(q)) * q**(a-1)
        Q, m = self.Q, self.m
        pf = list(factor(m))
        rng = random.Random(q*1000+a)
        while True:
            r = (rng.randrange(Q), rng.randrange(Q))
            if (r[0]*r[0]+r[1]*r[1]) % q == 0: continue
            t = ratio(r, Q)
            if all(gpow(t, m//l, Q) != (1, 0) for l in pf):
                self.g = t; break
        self.log = {}
        x = (1, 0)
        for k in range(m):
            self.log[x] = k; x = gmul(x, self.g, Q)
        assert len(self.log) == m
        self.iota = self.log[(0, 1)]
    def dlog_ratio(self, r):
        return self.log[ratio((r[0] % self.Q, r[1] % self.Q), self.Q)]
