"""Exact arithmetic in Z[i]/Q (Q odd prime power) used by the checks."""
from math import gcd
import random

def is_prime(n):
    if n < 2: return False
    small = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small:
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d//=2; s+=1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def primes_upto(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def mul(a, b, Q):
    return ((a[0]*b[0]-a[1]*b[1]) % Q, (a[0]*b[1]+a[1]*b[0]) % Q)

def conj(a, Q):
    return (a[0] % Q, (-a[1]) % Q)

def norm(a, Q):
    return (a[0]*a[0]+a[1]*a[1]) % Q

def inv(a, Q):
    n = norm(a, Q)
    ni = pow(n, -1, Q)
    c = conj(a, Q)
    return (c[0]*ni % Q, c[1]*ni % Q)

def pw(a, e, Q):
    r = (1, 0)
    if e < 0:
        a = inv(a, Q); e = -e
    while e:
        if e & 1: r = mul(r, a, Q)
        a = mul(a, a, Q); e >>= 1
    return r

def tq_order(q, a):
    chi = 1 if q % 4 == 1 else -1
    return (q - chi) * q**(a-1)

def v2(n):
    k = 0
    while n % 2 == 0: n//=2; k+=1
    return k

def sylow2_generator(q, a, rng):
    """element of T_Q of order exactly 2^k, k = v2(m_Q); deterministic given rng"""
    Q = q**a; m = tq_order(q, a); k = v2(m)
    while True:
        x = (rng.randrange(Q), rng.randrange(Q))
        if gcd(norm(x, Q), q) != 1: continue
        t = mul(x, inv(conj(x, Q), Q), Q)          # t in T_Q
        f = pw(t, m >> k, Q)                          # in 2-Sylow
        if pw(f, 1 << (k-1), Q) != (1, 0):            # order exactly 2^k
            return f, k

def gaussian_primes_norm_range(lo, hi):
    """one Gaussian prime x+iy (x>y>0) for each rational prime p=1 mod 4 in (lo,hi]"""
    out = []
    for p in primes_upto(hi):
        if p <= lo or p % 4 != 1: continue
        x = 1
        while True:
            y2 = p - x*x
            y = int(round(y2**.5))
            if y*y == y2:
                out.append((max(x,y), min(x,y))); break
            x += 1
    return out
