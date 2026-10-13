#!/usr/bin/env python3
"""Referee checks around Hypothesis H / Theorem C / Lemma C.3 (p3lower.md). Evidence + exact arithmetic.
1. Independent uniform classes on Z/m_q, window (M^(1-eta), M], eta = 1/2: Monte Carlo mean of F_k / E_k.
2. Exact e_k from the window primes; Lemma C.1 bounds.
3. Theorem C's k = floor((2 - eps/4) L / l): for which M is delta^k E_k >= 1 (eta=1/2, delta=1/2, eps=1/2)?
   Uses exact primes up to 10^7, then the Mertens value lambda = log 2 (approximation; labelled)."""
import math
import numpy as np
from sympy import primerange

rng = np.random.default_rng(99)

def mq(q): return q - 1 if q % 4 == 1 else q + 1

def esym(ps, kmax):
    e = [1.0] + [0.0] * kmax
    for p in ps:
        for k in range(kmax, 0, -1):
            e[k] += e[k - 1] * p
    return e

print("1. independent model, eta = 1/2, Monte Carlo (evidence)")
for M in [128, 512]:
    qs = [q for q in primerange(3, M + 1) if q > math.sqrt(M)]
    ps = [1.0 / mq(q) for q in qs]
    kmax = 6
    e = esym(ps, kmax)
    Ek = [M * (M - 1) / 2 * e[k] for k in range(kmax + 1)]
    acc = np.zeros(kmax + 1); T = 20
    for t in range(T):
        c = np.zeros((M, M), dtype=np.int64)
        for q in qs:
            kap = rng.integers(0, mq(q), size=M)
            c += (kap[:, None] == kap[None, :])
        iu = np.triu_indices(M, 1)
        cv = c[iu]
        for k in range(1, kmax + 1):
            acc[k] += sum(math.comb(int(v), k) for v in cv if v >= k)
    lam = sum(ps)
    print(f"  M={M}: lambda_Q={lam:.4f}; mean F_k/E_k, k=1..{kmax}: " +
          ", ".join(f"{acc[k]/T/Ek[k]:.3f}" for k in range(1, kmax + 1)))
    # Lemma C.1 bounds
    pmax = max(ps)
    for k in range(1, kmax + 1):
        lo = (lam - k * pmax) ** k / math.factorial(k); hi = lam ** k / math.factorial(k)
        assert lo <= e[k] <= hi + 1e-15, (k, lo, e[k], hi)
print("  Lemma C.1 bounds hold for the exact window e_k.")

print("3. Theorem C threshold: smallest M with delta^k E_k >= 1 for k = floor((2-eps/4)L/l)")
eta, delta, eps = 0.5, 0.5, 0.5
def logterm(M, lam, pmax):
    L = math.log(M); l = math.log(L)
    k = math.floor((2 - eps / 4) * L / l)
    lo = (lam - k * pmax)
    return k, math.log(M * (M - 1) / 2) + k * math.log(delta) + k * math.log(lo) - math.lgamma(k + 1)
for M in [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]:
    qs = [q for q in primerange(3, M + 1) if q > M ** (1 - eta)]
    lam = sum(1.0 / mq(q) for q in qs); pmax = 1.0 / mq(qs[0])
    k, lt = logterm(M, lam, pmax)
    print(f"  M=1e{round(math.log10(M))}: lambda_Q={lam:.4f}, k={k}, log(delta^k E_k) (Lemma C.1 lower) = {lt:.2f}")
lam = math.log(2)
for ex in [8, 10, 15, 20, 30, 40, 60, 100]:
    M = 10.0 ** ex
    L = math.log(M); l = math.log(L)
    k = math.floor((2 - eps / 4) * L / l)
    lt = 2 * L - math.log(2) + k * math.log(delta * lam) - math.lgamma(k + 1)
    print(f"  M=1e{ex} (Mertens lambda=log 2, approx): k={k}, log(delta^k E_k) ~ {lt:.2f}, guaranteed D_max >= 2(1-eta)kL = {2*(1-eta)*k*L:.1f} = {2*(1-eta)*k/L:.3f} L^2")
