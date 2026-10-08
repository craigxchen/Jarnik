#!/usr/bin/env python3
"""Evaluate the configuration-free finite function F*(M) of proof.md Theorem 5.2:

  F*(M) = binom(M,2) log2/2
        + sum_{l = 3 mod 4} log l * sum_{a>=1} E(M, (l+1) l^(a-1))
        + sum_{p = 1 mod 4} log p * G_p(M),
  G_p(M) = min( sum_{a>=1} E(M,(p-1)p^(a-1)),  Psi_p(M) ),
  Psi_p(M) = min over e>=1 and level distributions (n_0,n_e >= 1, no empty middle level) of
             sum_a sum_t E(n_t,(p-1)p^(a-1)) + (1/2) sum_tau (S_tau - M/2)^2.

Every configuration of M points satisfies F*(M) <= (M/4) log R' + binom(M,2) log C, with
R' <= R the normalised radius.  We print F*(M), the implied minimal radius for C = 1, and
the constant Kneed(M) = 3 log M - 8 F*(M)/M^2 (so that 3 M log M <= W + Kneed(M) M).
Integer arithmetic throughout (brackets are scaled by 8).
"""
from math import comb, log
import sys


def E(n, q):
    if n <= 0:
        return 0
    b, r = divmod(n, q)
    return q * comb(b, 2) + r * b


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def depth_sum(n, q0, p):
    tot = 0
    q = q0
    while q < n:
        tot += E(n, q)
        q *= p
    return tot


def psi8(M, p):
    """8 * Psi_p(M), exact integer."""
    w = [8 * depth_sum(n, p - 1, p) for n in range(M + 1)]
    INF = None
    f = [INF] * (M + 1)
    for S in range(1, M):
        f[S] = w[S] + (2 * S - M) ** 2
    for S in range(1, M):
        for S2 in range(S + 1, M):
            v = f[S] + w[S2 - S] + (2 * S2 - M) ** 2
            if v < f[S2]:
                f[S2] = v
    return min(f[S] + w[M - S] for S in range(1, M))


def Fstar(M, primes):
    tot = comb(M, 2) * log(2) / 2
    detail = {'inert': 0.0, 'split': 0.0, 'split_div_worse': 0}
    for q in primes:
        if q == 2:
            continue
        if q % 4 == 3:
            if q + 1 >= M:
                continue
            v = depth_sum(M, q + 1, q) * log(q)
            tot += v
            detail['inert'] += v
        else:
            if q - 1 >= M:
                continue
            free8 = 8 * depth_sum(M, q - 1, q)
            ps8 = psi8(M, q)
            g8 = min(free8, ps8)
            if ps8 < free8:
                detail['split_div_worse'] += 1
            v = g8 / 8 * log(q)
            tot += v
            detail['split'] += v
    return tot, detail


def main():
    Ms = [4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 160, 200]
    primes = primes_upto(max(Ms) + 2)
    print("    M        F*(M)    (3/8)M^2logM   inert part   split part  #p with Psi<free"
          "   Kneed(M)   logR_min(C=1)")
    for M in Ms:
        F, d = Fstar(M, primes)
        Kneed = 3 * log(M) - 8 * F / (M * M)
        logRmin = 4 * F / M
        print(f" {M:4d} {F:12.2f} {3/8*M*M*log(M):13.2f} {d['inert']:12.2f} {d['split']:11.2f}"
              f" {d['split_div_worse']:8d}        {Kneed:8.4f}   {logRmin:10.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
