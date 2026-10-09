"""Explicit size of the fake circles of Theorem 2 as a function of M (b = 8), from the proved bounds:
  * Theorem 1 margins g(2) = b-4, g(n) = min(2n+b-8, (b-3)n/2+b) (n >= 4);
  * |C_n| <= (2M)^n characters of l1-norm n;
  * residue design: L = M-th prime = 1 mod 4 above Q0 = 2M (class A1), or no residues (class A0, L = 1);
  * union bound  sum_n (2M)^n L^(n/2) (4*1.0368/C) e^(1/16) e^(-tau g(n)/4) * 2.01 <= 1/2;
  * prime window: P >= 8.4 r^2 log P (pigeonhole of the primes = 1 mod 4 in [P,2P] into windows of
    log-width 1/(4r), using pi(2P;4,1)-pi(P;4,1) >= P/(3 log P)), and p_j >= 21 r^2 (Theorem A.1(ii)).
  W <= r (log 2P + 1/(4r)).  Pure arithmetic in floating point (logs only); reports W/(M log M)."""
import math, sys

def primes_1mod4_above(Q0, count):
    hi = Q0 + 50 * count + 1000
    while True:
        s = bytearray([1]) * (hi + 1); s[0] = s[1] = 0
        for i in range(2, int(hi ** 0.5) + 1):
            if s[i]:
                s[i*i::i] = bytearray(len(s[i*i::i]))
        ps = [p for p in range(Q0 + 1, hi + 1) if s[p] and p % 4 == 1]
        if len(ps) >= count:
            return ps[:count]
        hi *= 2

def tau_union(M, b, C, lnL):
    g = lambda n: (b - 4) if n == 2 else min(2 * n + b - 8, (b - 3) * n / 2 + b)
    lo, hi = 0.0, 500.0
    for _ in range(100):
        tau = (lo + hi) / 2
        tot = 0.0
        # for n >= 4 (b = 8) g(n) = 2n, so the exponent is affine in n with slope
        # log(2M) + lnL/2 - tau/2; if the slope is >= 0 the series diverges, otherwise the terms decrease
        # geometrically and we stop once they are below e^-80.
        if math.log(2 * M) + lnL / 2 - tau / 2 >= 0:
            lo = tau
            continue
        n = 2
        while True:
            e = n * (math.log(2 * M) + lnL / 2) - tau * g(n) / 4 + 1 / 16
            tot += math.exp(min(e, 700)) * 4 * 1.0368 / C * 2.01
            if (n >= 8 and e < -80) or tot > 1:
                break
            n += 2
        if tot <= 0.5:
            hi = tau
        else:
            lo = tau
    return hi

def tau_window(r):
    P = 21.0 * r * r
    while P < 8.4 * r * r * math.log(P):
        P *= 1.01
    return math.log(P)

b = 8; C = 1.0
print(" M        r      tau_union(A1)  tau_window   W/(M log M) A1   W/(M log M) A0   tau/log M (A1)")
for M in [12, 20, 44, 104, 1004, 10004, 100004, 1000004]:
    r = b * (M - 1) + M
    L = primes_1mod4_above(2 * M, M)[-1] if M <= 100004 else None
    lnL = math.log(L) if L else math.log(2 * M + 3.2 * M * math.log(M))
    t1 = tau_union(M, b, C, lnL); t0 = tau_union(M, b, C, 0.0); tw = tau_window(r)
    W1 = r * (max(t1, tw) + math.log(2) + 1 / (4 * r)); W0 = r * (max(t0, tw) + math.log(2) + 1 / (4 * r))
    lm = M * math.log(M)
    print("%8d %8d   %8.2f      %8.2f      %8.2f          %8.2f          %6.3f" % (M, r, t1, tw, W1 / lm, W0 / lm, max(t1, tw) / math.log(M)))
