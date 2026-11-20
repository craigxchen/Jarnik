"""Finite-m thresholds for X.3: full-core column condition, 5/8 margin, deletion fraction.
Worst case: the t*m primes are the t*m consecutive primes starting at the smallest P with
pi(2P)-pi(P) >= t*m, columns get the largest primes in the worst place (we bound Q_v by (2P)^t
for the divisor and q >= P^{tm}; also report exact values with actual primes)."""
import math, sys
from itertools import combinations

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

PR = primes_upto(4_000_000)

def smallest_P(n):
    # smallest P such that [P,2P] has >= n primes
    import bisect
    lo, hi = 2, 2_000_000
    while lo < hi:
        mid = (lo + hi) // 2
        c = bisect.bisect_right(PR, 2 * mid) - bisect.bisect_left(PR, mid)
        if c >= n: hi = mid
        else: lo = mid + 1
    return lo

def is_prime(n):
    return n > 1 and all(n % p for p in range(2, int(n ** 0.5) + 1))

def data(m, t):
    import bisect
    P = smallest_P(t * m)
    i0 = bisect.bisect_left(PR, P)
    ps = PR[i0:i0 + t * m]
    logq = sum(math.log(p) for p in ps)
    return P, ps, logq

def full_core_ok(m, t, C):
    # need ((m+1)/2) * C * sqrt q < q/(Q_v Q_w) for worst v,w: Q_vQ_w <= product of 2t largest primes
    P, ps, logq = data(m, t)
    worst = sum(math.log(p) for p in ps[-2 * t:])
    lhs = math.log((m + 1) / 2) + math.log(C) + logq / 2 + worst
    # triangle step: g_min >= product of t(m-3)/4 smallest primes ... need g > C^2
    g_ok = (m - 3) // 4 >= 1 and t * ((m - 3) // 4) * math.log(P) > 2 * math.log(C)
    return lhs < logq, g_ok, P

def deletion_kmax(m, t, C, sharp=False):
    P, ps, logq = data(m, t)
    best = -1
    for k in range(0, m - 4):
        r = k + 1
        worst = sum(math.log(p) for p in ps[-t * (r + 1):])  # largest t(r+1) primes
        if sharp:
            # |a^t y| <= (||a||_1 / 2) L, ||a||_1 <= m (r+1) r^{r/2}
            lhs = math.log(m * (r + 1) / 2) + (r / 2) * math.log(r) + math.log(C) + logq / 2 + worst
        else:
            lhs = math.log(m * (r + 1)) + (r / 2) * math.log(r) + math.log(C) + logq / 2 + worst
        if lhs < logq:
            best = k
        else:
            break
    return best, P

def margin(m, t):
    # min over triples of log g / log q (equal t per column); worst case: smallest primes in constant cols
    P, ps, logq = data(m, t)
    c = (m - 3) // 4
    gmin = sum(math.log(p) for p in ps[:t * c])
    return gmin / logq

if __name__ == "__main__":
    ms = [m for m in range(7, 20000) if m % 4 == 3 and is_prime(m)]
    for t in (1, 2):
        for C in (1, 10, 1000):
            m0 = next((m for m in ms if all(full_core_ok(mm, t, C)[0] and full_core_ok(mm, t, C)[1]
                                         for mm in ms if mm >= m and mm <= m + 200)), None)
            print(f"t={t} C={C}: full-core column+triangle condition holds from m={m0} (checked to m+200)")
    for t in (1, 2, 3):
        for m in (43, 103, 503, 1999, 4999, 9999 if is_prime(9999) else 9973 if 9973 % 4 == 3 else 10007, 19991 if False else 19963):
            if not (is_prime(m) and m % 4 == 3):
                continue
            k, P = deletion_kmax(m, t, 1)
            ks, _ = deletion_kmax(m, t, 1, sharp=True)
            print(f"t={t} m={m} P={P}: max deletions (her ineq, C=1) k={k}, k/m={k/m:.3f}; "
                  f"retained fraction {(m-k)/m:.3f}; asymptotic k/m = {t/(2*t+1):.3f}; margin logg/logq={margin(m,t):.4f}")
            sys.stdout.flush()
