#!/usr/bin/env python3
"""Families near alpha = 1/4.
(a) n = prod_{i<=4}(1+m_i q), arbitrary (composite allowed) factors: count divisors of n in
    each invertible class mod q exactly; report best counts with alpha >= theta.
(b) the q=1 interval version: n = x(x+a)(x+b)(x+a+b) has the 6 pair products in an interval
    [A,B] with (B-A)/sqrt(AB) ~ c n^{-1/4}."""
import random
from math import log, gcd, isqrt
from sympy import divisors, factorint
rng = random.Random(11)
thetas = [0.20, 0.22, 0.23, 0.24, 0.245, 0.25]
best = {t: (0, None) for t in thetas}
for trial in range(4000):
    q = rng.randint(30, 3000)
    k = rng.choice([4, 4, 4, 5])
    ms = [rng.randint(1, 40) for _ in range(k)]
    fac = [1 + m * q for m in ms]
    n = 1
    for f in fac: n *= f
    al = log(q) / log(n)
    if al < thetas[0]: continue
    cc = {}
    for d in divisors(n):
        if gcd(d, q) == 1:
            cc[d % q] = cc.get(d % q, 0) + 1
    r, c = max(cc.items(), key=lambda kv: kv[1])
    for t in thetas:
        if al >= t and c > best[t][0]:
            best[t] = (c, (q, ms, n, r, al))
for t in thetas:
    print("theta", t, best[t])

print("(b) q=1 interval family")
for x in [10**3, 10**6, 10**12, 10**24]:
    a, b = 1, 2
    fs = [x, x + a, x + b, x + a + b]
    n = fs[0] * fs[1] * fs[2] * fs[3]
    prods = sorted({fs[i] * fs[j] for i in range(4) for j in range(i + 1, 4)})
    assert all(n % d == 0 for d in prods)
    A, B = prods[0], prods[-1]
    ratio = (B - A) / (A * B) ** 0.5
    print(f"x=1e{round(log(x,10))}: {len(prods)} divisors, (B-A)/sqrt(AB) * n^(1/4) = {ratio * n**0.25:.4f}")
