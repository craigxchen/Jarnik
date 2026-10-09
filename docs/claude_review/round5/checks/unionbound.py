"""Evaluates the union bound of Theorem 5.1(b) of lower.md (Lemma S route, any capacity-two
assignment) for explicit parameters.
ILLUSTRATION of the size of M_0 only: the parameter choices P = b^3 M^2 (log M)^2 and
Y = 3 M log M presume enough primes = 1 mod 4 in [P, 2P] and in (2M, Y] (true for large M by the
prime number theorem for arithmetic progressions; not certified here).  The proof of Theorem 5.1 is
asymptotic and does not use this script.

Bound:  sum over characters (mod +-) of (4 pi Gamma / C) (1 + n) m_c exp(-s_c / 2), where
  pairs        s >= (b/2 - 2) tau - 2 log 2 - 1/4,        log m <= 4 log Y - log 4
  columns      s >= (b/2 + M - 4) tau - 2 log 2 - 1/2,    log m <= L_X
  half-sums    s >= (b/2 + M/2 - 8) tau - 4 log 2 - 1/2,  log m <= L_X
  others       s >= (b tau/2) max(theta M, n - M) - (tau + log 2) n - 1/2,
               log m <= min(2 n log Y, L_X),
with L_X = log lcm(1..X) over odd primes, X = 2M, theta = 1/16, and #{c : ||c||_1 = n} computed
by C(n+M-1,M-1) 2^min(n,M) (no zero-sum or +- reduction: an over-count)."""
import math


def log_lcm_odd(X):
    s = 0.0
    sieve = bytearray([1]) * (X + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(X ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    for q in range(3, X + 1):
        if sieve[q]:
            a = 1
            while q ** (a + 1) <= X:
                a += 1
            s += a * math.log(q)
    return s


def log_count(M, n):
    # log of the upper bound  #{c in Z^M : ||c||_1 = n} <= C(n+M-1, M-1) 2^min(n,M)
    return (math.lgamma(n + M) - math.lgamma(M) - math.lgamma(n + 1)
            + min(n, M) * math.log(2))


def logsumexp(xs):
    mx = max(xs)
    return mx + math.log(sum(math.exp(x - mx) for x in xs))


def total(M, b, C=1.0, Gamma=2.0, theta=1 / 16):
    L = math.log(M)
    P = b ** 3 * M ** 2 * L ** 2
    tau = math.log(P)
    Y = 3 * M * L
    LX = log_lcm_odd(2 * M)
    pre = math.log(4 * math.pi * Gamma / C)
    parts = {}
    s = (b / 2 - 2) * tau - 2 * math.log(2) - 0.25
    parts["pairs"] = math.log(M * M / 2) + pre + math.log(3) + 4 * math.log(Y) - math.log(4) - s / 2
    s = (b / 2 + M - 4) * tau - 2 * math.log(2) - 0.5
    parts["columns"] = math.log(M) + pre + math.log(1 + M) + LX - s / 2
    s = (b / 2 + M / 2 - 8) * tau - 4 * math.log(2) - 0.5
    parts["half-sums"] = math.log(M * M) + pre + math.log(1 + M / 2) + LX - s / 2
    xs = []
    nmax = 40 * M
    for n in range(4, nmax + 1, 2):
        s = (b * tau / 2) * max(theta * M, n - M) - (tau + math.log(2)) * n - 0.5
        lm = min(2 * n * math.log(Y), LX)
        xs.append(log_count(M, n) + pre + math.log(1 + n) + lm - s / 2)
    parts["others"] = logsumexp(xs)
    tot = logsumexp(list(parts.values()))
    W_over = b * (tau + math.log(2)) * (M - 1) / (M * L)
    return tot, parts, W_over


if __name__ == "__main__":
    for b in (35, 40, 48, 65):
        for M in (44, 200, 1000, 5000):
            tot, parts, W_over = total(M, b)
            print("b=%3d M=%6d  log(total bound)=%10.2f  [pairs %.1f, columns %.1f, half %.1f, others %.1f]"
                  "  W/(M log M) = %.1f" % (b, M, tot, parts["pairs"], parts["columns"],
                                            parts["half-sums"], parts["others"], W_over))
