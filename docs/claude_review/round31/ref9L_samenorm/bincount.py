#!/usr/bin/env python3
"""Count split primes (p = 1 mod 4) in one bin [P, P*exp(l^{-3/4})], P = l^1.99, by a segmented sieve.
Compare with the PNT prediction for the average bin: (#split primes in [P,2P]) / (#bins),
#bins = ceil(log 2 * l^{3/4}). Also report K/s, J = log(K/s), log(2P)/J with s = l-1."""
import sys, math
import numpy as np
from sympy import isprime, nextprime

def base_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def count_split(lo, hi, chunk=2*10**7):
    bp = base_primes(int(math.isqrt(hi)) + 1)
    total = 0
    a = lo
    while a < hi:
        b = min(hi, a + chunk)
        seg = np.ones(b - a, dtype=bool)
        for p in bp:
            p = int(p)
            if p * p >= b: break
            st = max(p * p, ((a + p - 1) // p) * p)
            seg[st - a::p] = False
        idx = np.nonzero(seg)[0] + a
        total += int(np.count_nonzero(idx % 4 == 1))
        a = b
    return total

def li(x):
    from mpmath import li as mli
    return float(mli(x))

for l0 in [10**4, 10**5, 10**6, 10**7]:
    l = l0
    while not (isprime(l) and l % 4 == 3): l += 1
    P = int(round(l ** 1.99))
    sig = l ** -0.75
    hi = int(P * math.exp(sig))
    K = count_split(P, hi)
    nb = math.ceil(math.log(2) / sig)
    avg_pred = (li(2*P) - li(P)) / 2 / nb
    asym = l ** 1.24 / (2 * 1.99 * math.log(2) * math.log(l))
    s = l - 1
    J = math.log(K / s) if K > 0 else float('nan')
    print(f"l={l}: P={P:.3e}, bin length={hi-P}, K(first bin)={K}, PNT avg-bin={avg_pred:.0f}, "
          f"l^1.24/(2*1.99*log2*log l)={asym:.0f}, s={s}, K/s={K/s:.3f}, J={J:.3f}, "
          f"log(2P)/J={math.log(2*P)/J if K> s else float('nan'):.2f}", flush=True)
# asymptotic thresholds
for name, J0 in [("J>=40 (squarefree IX.0)", 40.0), ("J>=400E^2 log(E(E+1)^3), E=2 (ref9L_powers Thm G)", 1600*math.log(54))]:
    # solve 0.24 log l - log(2*1.99*log2*log l) = J0
    L = 10.0
    for _ in range(200):
        L = (J0 + math.log(2*1.99*math.log(2)*L)) / 0.24
    print(f"{name}: needs log l >= about {L:.0f}")
