#!/usr/bin/env python3
"""Numerical constants for proof.md, Lemma 7.4 / Theorem 7.5, with rigorous tail bounds, plus sanity
checks of the effective inequality against the exact finite quantities.

  sigma = sum_p log p/(p(p-1))
  a3    = sum_{l = 3 mod 4} log l/(l(l+1))
  a1    = sum_{p = 1 mod 4} log p/(p(p-1)^2)
Each is computed as a partial sum over p <= P plus an explicit upper bound for the tail
(f(t) = log t/(t(t-1)) is decreasing for t >= 2, so sum_{n>P} f(n) <= int_P^oo f
<= (P/(P-1)) (log P + 1)/P).  Floating-point error is far below the 1e-9 safety margin.
"""
from math import log, pi, e, sqrt, comb
import sys

P = 10 ** 7


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def main():
    primes = sieve(P)
    tail = (P / (P - 1)) * (log(P) + 1) / P
    margin = 1e-9
    sigma = sum(log(p) / (p * (p - 1)) for p in primes)
    a3 = sum(log(p) / (p * (p + 1)) for p in primes if p % 4 == 3)
    a1 = sum(log(p) / (p * (p - 1) ** 2) for p in primes if p % 4 == 1)
    sigma_up = sigma + tail + margin
    a3_up = a3 + tail + margin
    a1_up = a1 + tail / (P - 1) + margin
    print(f"sigma  in [{sigma:.10f}, {sigma_up:.10f}]")
    print(f"a3     in [{a3:.10f}, {a3_up:.10f}]")
    print(f"a1     in [{a1:.10f}, {a1_up:.10f}]")
    beta = 1 + sigma_up + log(2) / 2
    xi = (4 / pi) * (log(4) + 2 / e)
    c0 = 2 * a3_up + a1_up
    K0 = 3 * beta + xi + 2 * c0 + 6 * log(2)
    print(f"beta = 1 + sigma + log2/2      <= {beta:.6f}")
    print(f"xi   = (4/pi)(log 4 + 2/e)     <= {xi:.6f}")
    print(f"c0   = 2 a3 + a1               <= {c0:.6f}")
    print(f"6 log 2                         = {6*log(2):.6f}")
    print(f"K0 = 3 beta + xi + 2 c0 + 6 log 2 <= {K0:.6f}")

    # ---------------------------------------------------------------
    # sanity check 1: depth-one analytic lower bound D(M) of (6.2), C = 1:
    #   D(M) = sum_{l<=M, l=3(4)} (M^2/(2(l+1)) - M/2) log l
    #        + sum_{p<=M, p=1(4)} (c_p M^2/4 - M/2) log p + binom(M,2) log2/2,
    # the finite inequality gives (M/8) W >= D(M), i.e. 3M log M - W <= M*Kneed(M),
    #   Kneed(M) = 3 log M - 8 D(M)/M^2.   Must be <= K0 for every M (Theorem 7.5).
    # computed incrementally for all M up to MMAX.
    MMAX = 2 * 10 ** 6
    ps = [p for p in primes if p <= MMAX]
    import bisect
    A2 = 0.0   # sum log l/(l+1) over inert l <= M
    B2 = 0.0   # sum c_p log p over split p <= M
    TH = 0.0   # sum log q over odd q <= M
    idx = 0
    worst = (-1e9, None)
    rows = []
    checkpoints = {10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 2 * 10 ** 6}
    for M in range(2, MMAX + 1):
        while idx < len(ps) and ps[idx] <= M:
            q = ps[idx]
            if q % 4 == 3:
                A2 += log(q) / (q + 1)
                TH += log(q)
            elif q % 4 == 1:
                lam = 1 / (q - 1)
                cq = 2 * lam / (sqrt(1 + 4 * lam) + 1)  # = (sqrt(1+4lam)-1)/2, stable
                B2 += cq * log(q)
                TH += log(q)
            idx += 1
        D = M * M / 2 * A2 + M * M / 4 * B2 - M / 2 * TH + comb(M, 2) * log(2) / 2
        Kneed = 3 * log(M) - 8 * D / (M * M)
        if Kneed > worst[0]:
            worst = (Kneed, M)
        if M in checkpoints:
            rows.append((M, Kneed, 2 * A2 + B2 - 1.5 * log(M)))
    print("\nDepth-one analytic bound D(M):  Kneed(M) = 3 log M - 8 D(M)/M^2  (C = 1)")
    print("      M      Kneed(M)    2A(M)+B(M)-(3/2)log M")
    for r in rows:
        print(f" {r[0]:8d}  {r[1]:10.4f}  {r[2]:10.4f}")
    print(f"max_{{2<=M<={MMAX}}} Kneed(M) = {worst[0]:.4f} at M = {worst[1]}  (proved bound K0 = {K0:.4f})")
    ok = worst[0] <= K0
    print("Kneed(M) <= K0 on the whole range:", ok)

    # ---------------------------------------------------------------
    # sanity check 2 (numerical illustration of Lemmas 7.1-7.3 on n <= 2e6; they are proved):
    #   O(n) = sum_{3<=p<=n} log p/p >= log n - beta + 1/n,  X(n) = sum chi(p) log p/p <= xi,
    #   theta(n) < n log 4.
    Osum = 0.0; Xsum = 0.0; th = 0.0; idx = 0
    minO = 1e9; maxX = -1e9; maxth = 0.0; maxpsi_ratio = 0.0
    for n in range(2, MMAX + 1):
        while idx < len(ps) and ps[idx] <= n:
            q = ps[idx]; th += log(q)
            if q > 2:
                Osum += log(q) / q
                Xsum += (1 if q % 4 == 1 else -1) * log(q) / q
            idx += 1
        minO = min(minO, Osum - (log(n) - beta + 1 / n))
        maxX = max(maxX, Xsum)
        maxth = max(maxth, th / n)
    print(f"\nmin_n [O(n) - (log n - beta + 1/n)] = {minO:.4f} (>= 0 required)")
    print(f"max_n X(n) = {maxX:.4f}  (proved bound xi = {xi:.4f})")
    print(f"max_n theta(n)/n = {maxth:.4f}  (proved bound log 4 = {log(4):.4f})")
    ok = ok and minO >= 0 and maxX <= xi and maxth < log(4)

    # ---------------------------------------------------------------
    # explicit threshold of Corollary 8.1:  ell = loglog R must satisfy
    #   (i)  (log ell + log(3/2) + K/3)/ell <= 3 eps/(2+3 eps),   (ii) (2/3) e^ell/ell >= e^(K/3)
    from math import exp
    print("\nThreshold ell_0 = loglog R_0 of Corollary 8.1 (K = K_C):")
    for K in (14.0, 14.0 + 4 * log(10)):
        for eps in (1 / 3, 0.1, 0.01):
            ell = 1.0
            while True:
                c1 = (log(ell) + log(1.5) + K / 3) / ell <= 3 * eps / (2 + 3 * eps)
                c2 = log(2 / 3) + ell - log(ell) >= K / 3
                if c1 and c2:
                    break
                ell += 0.01
            print(f"  K={K:7.3f} eps={eps:6.3f}:  ell_0 <= {ell:9.2f}   (R_0 = exp(exp(ell_0)))")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
