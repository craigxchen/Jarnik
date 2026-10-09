"""Referee (analytic lens): Theorem 6.4 (P3-all) step 4, re-derived and tested.

(i)  per-prime moment bound. If Pr[a_q >= a] <= t_a := min(1, 24 kk q^-a) (a <= A_q), the largest
     possible E[q^(theta a_q)] is 1 + sum_a t_a (q^(theta a) - q^(theta(a-1))) (tail-dominated laws).
     Check against sharp.md:  q > 24 kk: <= 1 + 58 kk q^(theta-1);   q <= 24 kk: <= 1 + (log_3(24X)+3.4) e^(2 beta),
     with X^theta = e^beta; grid over q, kk, theta in (0, 1/2], X.
(ii) the summed bound  sum_(q <= X) log E_max(q)  <=  kk [20 beta + 10 loglog X + 10 + 58 e^beta (loglog X + 1)]
     with exact prime sums for X up to 10^7.
(iii) the count #{v in Z^r : ||v||_1 = kk} <= (2r)^kk (exact for small r, kk).
(iv) parameters: e^beta = log(4er)/(58 (loglog X + 1) loglog M), theta = beta/log X, tau_* ; where beta > 0
     starts, and tau_*/(log^2 M/loglog M) against the claimed limit 2A.
"""
import math
import numpy as np

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]

def Emax(q, kk, theta, X):
    Aq = int(math.floor(math.log(X) / math.log(q) + 1e-12))
    while q ** (Aq + 1) <= X: Aq += 1
    while Aq > 0 and q ** Aq > X: Aq -= 1
    E = 1.0
    for a in range(1, Aq + 1):
        t = min(1.0, 24 * kk * q ** (-a))
        E += t * (q ** (theta * a) - q ** (theta * (a - 1)))
    return E

worst1 = worst2 = 0.0
for X in (10 ** 3, 10 ** 6, 10 ** 12, 10 ** 30):
    for theta in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5):
        beta = theta * math.log(X)
        for kk in (1, 2, 5, 20, 100, 10 ** 4):
            for q in (3, 5, 7, 11, 31, 101, 499, 1009, 10007, 10 ** 5 + 3, 10 ** 6 + 3):
                if q > X:
                    continue
                E = Emax(q, kk, theta, X)
                if q > 24 * kk:
                    worst1 = max(worst1, (E - 1) / (58 * kk * q ** (theta - 1)))
                else:
                    worst2 = max(worst2, (E - 1) / ((math.log(24 * X) / math.log(3) + 3.4) * math.exp(2 * beta)))
print(f"(i)  max (E_max - 1)/(58 kk q^(theta-1)) over q > 24kk: {worst1:.3f}  (<= 1 required)")
print(f"     max (E_max - 1)/((log_3(24X)+3.4) e^(2 beta)) over q <= 24kk: {worst2:.3f}  (<= 1 required; trivial bound X^theta = e^beta also works)")

worst = 0.0
for X in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7):
    P = [int(q) for q in primes_upto(X) if q > 2]
    llX = math.log(math.log(X))
    for beta in (0.25, 0.5, 1.0, 2.0, 3.0):
        theta = beta / math.log(X)
        if theta > 0.5:
            continue
        for kk in (1, 3, 10, 100, 10 ** 4):
            s = sum(math.log(Emax(q, kk, theta, X)) for q in P)
            claim = kk * (20 * beta + 10 * llX + 10 + 58 * math.exp(beta) * (llX + 1))
            worst = max(worst, s / claim)
print(f"(ii) max over grid of sum_q log E_max(q) / claimed kk-bracket: {worst:.3f}  (<= 1 required)")

okc = True
from itertools import product
for r in (1, 2, 3, 4):
    for kk in (1, 2, 3, 4, 5):
        cnt = sum(1 for v in product(range(-kk, kk + 1), repeat=r) if sum(abs(x) for x in v) == kk)
        okc &= cnt <= (2 * r) ** kk
print(f"(iii) #{{||v||_1 = kk}} <= (2r)^kk for r <= 4, kk <= 5: {okc}")

print("(iv) parameters (b = 33, r = 33 M):")
for A in (1, 2):
    first = None
    for L in np.arange(10, 20000, 1.0):

        X = math.log(3) + A * L                      # log X
        eb = (math.log(4 * math.e) + math.log(33) + L) / (58 * (math.log(X) + 1) * math.log(L))
        if eb > 1 and first is None:
            first = L
    print(f"   A = {A}: beta > 0 first at log M = {first:.0f}")
    for L in (1e4, 1e5, 1e6, 1e9, 1e15):
        logr = math.log(33) + L
        logX = math.log(3) + A * L
        ll = math.log(L); llX = math.log(logX)
        eb = (math.log(4 * math.e) + logr) / (58 * (llX + 1) * ll)
        beta = math.log(eb)
        tau = (2 * logX / beta) * (math.log(4 * math.e) + logr + 14 + 20 * beta + 10 * llX + 58 * eb * (llX + 1))
        print(f"      log M = {L:9.0e}: beta = {beta:7.3f}, theta = {beta / logX:.2e}, tau_*/(log^2 M/loglog M) = {tau / (L * L / ll):7.3f}  (limit {2 * A})")
