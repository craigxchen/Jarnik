"""Referee check of the constants in the proof of Theorem A.3 (outside.md, Section 2.2).

(a) Lemma A.6 hypothesis n > k^2 m (k log 2/eta + 1), k = 6M-4, m = M-1, eta = 1/M, plus 4(M-1) primes
    for F and S: is P = 400 M^5 log M enough?  Exact count of primes = 1 mod 4 in [P, 2P] for M = 4, 8
    (segmented numpy sieve), PNT estimate P/(2 log P) for larger M.
(b) The smallest P (power-of-1.1 grid) for which the PNT estimate meets the requirement, as a multiple of
    M^5 log M.
(c) The union-bound total of step 3, 1.04 (2.55 + 8/C) (3 e^2.08 P^-0.43)^M sum_{n>=2} n (2/sqrt P)^n / 2,
    for C = 1/2, 1e-3 at the corrected P.
(d) Q_M = lcm(1..2M) <= e^(2.08 M) for M <= 2000 (exact integers).
(e) arcsin y <= 1.0367 y on (0, 5^-1/2], and (8/pi)*1.0367/2 <= 1.32.
"""
from math import log, exp, sqrt, asin, pi, lcm
import numpy as np


def count_1mod4(lo, hi):
    """number of primes p = 1 mod 4 with lo <= p <= hi (sieve)."""
    n = hi + 1
    s = np.ones(n, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    idx = np.nonzero(s[lo:])[0] + lo
    return int(np.count_nonzero(idx % 4 == 1))


def need(M):
    k, m, eta = 6 * M - 4, M - 1, 1.0 / M
    return k * k * m * (k * log(2) / eta + 1) + 4 * (M - 1)


print("(a) Lemma A.6 requirement vs primes = 1 mod 4 in [P,2P] at the text's P = 400 M^5 log M")
for M in (4, 8):
    P = int(400 * M ** 5 * log(M))
    c = count_1mod4(P, 2 * P)
    print("   M=%3d  P=%.3e  required > %.3e  available (exact sieve) %.3e  -> %s"
          % (M, P, need(M), c, "OK" if c > need(M) else "INSUFFICIENT"))
for M in (16, 32, 64, 128, 1024):
    P = 400 * M ** 5 * log(M)
    est = P / (2 * log(P))
    print("   M=%4d  P=%.3e  required > %.3e  available (PNT est.) %.3e  ratio %.2f"
          % (M, P, need(M), est, est / need(M)))

print("(b) smallest grid P meeting the requirement (PNT estimate), in units of M^5 log M")
for M in (8, 16, 32, 64, 128, 1024, 10 ** 5):
    P = 1e3
    while P / (2 * log(P)) <= need(M):
        P *= 1.1
    print("   M=%6d  P/(M^5 log M) = %.0f" % (M, P / (M ** 5 * log(M))))
M = 8
P = int(1700 * M ** 5 * log(M))
c = count_1mod4(P, 2 * P)
print("   exact sieve M=8, P=1700 M^5 log M=%.3e: available %d vs required %.0f -> %s"
      % (P, c, need(8), "OK" if c > need(8) else "INSUFFICIENT"))

print("(c) union bound of step 3 at P = 1700 M^5 log M")
for C in (0.5, 1e-3):
    for M in (4, 8, 12, 16, 32):
        P = 1700 * M ** 5 * log(M)
        x = 2 / sqrt(P)
        tail = sum(n * x ** n for n in range(2, 200)) / 2
        val = 1.04 * (2.55 + 8 / C) * (3 * exp(2.08) * P ** -0.43) ** M * tail
        print("   C=%g M=%3d tau=log P=%.1f  bound=%.3e  %s" % (C, M, log(P), val, "<1" if val < 1 else ">=1"))

print("(d) lcm(1..2M) <= e^(2.08 M):", end=" ")
L = 1
ok = True
worst = 0
for M in range(1, 2001):
    L = lcm(L, 2 * M - 1, 2 * M)
    ratio = log(L) / M
    worst = max(worst, ratio)
    ok &= ratio <= 2.08
print("OK" if ok else "FAIL", " max log(Q_M)/M = %.4f" % worst)

ys = np.linspace(1e-9, 5 ** -0.5, 200001)
r = max(asin(y) / y for y in ys)
print("(e) max arcsin(y)/y on (0,5^-1/2] = %.5f ; (8/pi)*%.4f/2 = %.4f" % (r, r, 8 / pi * r / 2))

print("(f) referee extension: step-3 union bound with per-character strengthening Q_M * M^(2n), n=||c||_1,")
print("    i.e. sum_n 3^M 2^(n-1) 1.04 Q_M M^(2n) P^-(n+0.86M)/2 (2.55 n + 8/C), at P = 2000 M^5 log M")
for C in (0.5, 1e-3):
    for M in (4, 8, 16, 32, 64):
        P = 2000 * M ** 5 * log(M)
        x = 2 * M * M / sqrt(P)
        tail = sum(n * x ** n for n in range(2, 400)) / 2
        val = 1.04 * (2.55 + 8 / C) * (3 * exp(2.08) * P ** -0.43) ** M * tail
        print("   C=%g M=%3d  2M^2/sqrt(P)=%.3f  bound=%.3e  %s" % (C, M, x, val, "<1" if val < 1 else ">=1"))
