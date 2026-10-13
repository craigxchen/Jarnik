"""Referee ILLUSTRATION (float): Theorem B / Lemma 3.2 are proved for PRIME-supported inputs only.
The standard Hecke large sieve is stated over all ideals (all Gaussian integers mod units).  A fake has to be
extended to composites; the natural extension is multiplicative: every n divisible by varpi_p (p in S) has its
angle shifted by v_{varpi_p}(n) (phi_p - theta_p) - v_{conj}(n)(...).  Lemma 3.2's Cauchy-Schwarz bound for the
change is ~ (2K+1) * 4 * #{n : S | n} ||a||^2 ~ K X s/P_0, which is NOT o(K X^(1/2)) once X >> (P_0/s)^2.
Here we compute the optimal constant Lambda(X,K) (largest Gram eigenvalue) for all Gaussian integers of norm <= X,
first quadrant (one per associate class), genuine versus multiplicatively grafted at one or two primes, to see
whether the true change is small (it is, here) although the Lemma 3.2 bound is not.
"""
import math
import numpy as np
from sympy import isprime, factorint

X = 3000
pts = [(x, y) for x in range(1, math.isqrt(X) + 1) for y in range(0, math.isqrt(X) + 1) if 0 < x*x + y*y <= X]
ang = np.array([math.atan2(y, x) for x, y in pts])

def gauss_val(x, y, a, b):
    """valuation of x+iy at the Gaussian prime a+bi (p = a^2+b^2 split)."""
    v = 0; p = a*a + b*b
    while True:
        # (x+iy)/(a+bi) = (x+iy)(a-bi)/p
        re, im = x*a + y*b, y*a - x*b
        if re % p or im % p: return v
        x, y = re // p, im // p; v += 1

def grafted(primes_shift):
    out = ang.copy()
    for (a, b, d) in primes_shift:
        for i, (x, y) in enumerate(pts):
            v1 = gauss_val(x, y, a, b); v2 = gauss_val(x, y, b, a)   # varpi = a+bi and i*conj = b+ai
            out[i] += (v1 - v2) * d
    return out

def lam(angles, K):
    psi = 4 * angles
    d = psi[:, None] - psi[None, :]
    s = np.sin(d / 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        D = np.sin((K + 0.5) * d) / s
    D[np.abs(s) < 1e-13] = 2 * K + 1
    return np.linalg.eigvalsh(D).max()

rng = np.random.default_rng(1)
shifts = [(10, 1, rng.uniform(0, math.pi / 2)), (14, 1, rng.uniform(0, math.pi / 2))]   # p = 101, 197
g1 = grafted(shifts[:1]); g2 = grafted(shifts)
cnt = sum(1 for (x, y) in pts if any(gauss_val(x, y, a, b) or gauss_val(x, y, b, a) for a, b, _ in shifts))
print("all Gaussian integers (mod units) of norm <= %d: %d points; S-divisible (S = {101, 197}): %d" % (X, len(pts), cnt))
for K in (10, 30, 100, 300):
    l0, l1, l2 = lam(ang, K), lam(g1, K), lam(g2, K)
    crude = (2*K + 1) * 4 * cnt
    print("K=%4d: Lambda genuine %.0f, grafted@101 %.0f, grafted@101,197 %.0f; X + K X^(1/2) = %.0f; "
          "Lemma-3.2-type bound on the change (2K+1)*4*#S-div = %.0f" % (K, l0, l1, l2, X + K * math.sqrt(X), crude))
