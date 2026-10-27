"""Count S-integers (S = Gaussian primes over split p <= P, with the 4 units) in the annulus R/2 < |w| <= R,
and the mean number mu per square of side floor(sqrt(R)/2). Chebyshev: #hyperplanes needed >= ~mu/5.
Run: python3 -I count_mu.py"""
import math, sys, time
def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
def count(P, R):
    X = R * R; lo = X // 4
    ps = [p for p in primes_upto(P) if p % 4 == 1]
    gp = sorted([p for p in ps for _ in (0, 1)])  # each split p gives 2 Gaussian primes (as norms)
    cnt = 0
    sys.setrecursionlimit(10000)
    stack = [(0, 1)]
    while stack:
        k, n = stack.pop()
        if n > lo: cnt += 1
        for idx in range(k, len(gp)):
            q = n * gp[idx]
            if q > X: break
            stack.append((idx, q))
    U = 4 * cnt
    side = math.isqrt(R) // 2
    ncell = math.pi * (X - lo) / side ** 2
    return len(ps), U, U / ncell
for P, R in [(300, 3*10**4), (300, 10**5), (600, 10**5), (1000, 10**5)]:
    t = time.time(); s, U, mu = count(P, R)
    print(f"P={P} s={s} R={R}: U={U}  mean mu per square={mu:.1f}  mu/s={mu/s:.2f}  ({time.time()-t:.0f}s)", flush=True)
