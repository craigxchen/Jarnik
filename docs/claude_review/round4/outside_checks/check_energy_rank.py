"""Exact checks for Prop B.2 (four-term characters are subcritical on average) and the
prime-count step of Prop A.4 in outside.md.

B.2: for a 0/1 column with S ones among M rows, summing over ALL ordered quadruples (a,b,c,d)
     in [M]^4 (repetitions allowed):
        sum |x_a + x_b - x_c - x_d| = 4 S (M-S) (M^2 - S M + S^2)      (exact identity)
     while the pair sum over ordered pairs is  sum |x_a - x_c| * M^2 = 2 S (M-S) M^2.
     Hence quadruple average / pair average = 2 (1 - s + s^2) >= 3/2, s = S/M.
     For multi-level columns (values in {0..e}) we check the weaker exact inequality
        E|D + D'| >= E|D|   (D, D' iid copies of x_a - x_c)
     and record the minimal ratio found.
A.4: sum_{j<=k} log p_j >= k log(4k/e) for the first k primes p_j = 1 mod 4 (proved in the text
     from p_j >= 4j+1 and k! >= (k/e)^k); checked for k <= 10^5.
"""
import random
from fractions import Fraction as Fr
from math import log, e
from itertools import product

random.seed(7)


def quad_sum_binary(M, S):
    xs = [1] * S + [0] * (M - S)
    # use counts: sum over (i,j,k,l) of |xi+xj-xk-xl| via distribution of xi+xj
    cnt = {0: (M - S) ** 2, 1: 2 * S * (M - S), 2: S * S}
    return sum(cnt[u] * cnt[v] * abs(u - v) for u in cnt for v in cnt)


def check_binary():
    for M in range(1, 40):
        for S in range(0, M + 1):
            assert quad_sum_binary(M, S) == 4 * S * (M - S) * (M * M - S * M + S * S)
    print("B.2 binary identity: exact for all 1 <= M < 40, 0 <= S <= M: OK")


def check_multilevel(trials=3000):
    worst = None
    for _ in range(trials):
        e_ = random.randint(1, 5)
        M = random.randint(2, 12)
        xs = [random.randint(0, e_) for _ in range(M)]
        # distribution of D = xa - xc over ordered pairs
        D = {}
        for a in xs:
            for c in xs:
                D[a - c] = D.get(a - c, 0) + 1
        ED = sum(abs(d) * m for d, m in D.items())          # times M^2
        EDD = sum(abs(d1 + d2) * m1 * m2 for d1, m1 in D.items() for d2, m2 in D.items())  # times M^4
        if ED == 0:
            continue
        ratio = Fr(EDD, ED * M * M)
        assert ratio >= 1, (xs, ratio)
        worst = ratio if worst is None else min(worst, ratio)
    print("B.2 multilevel: E|D+D'| >= E|D| on %d random columns; min ratio = %s = %.4f"
          % (trials, worst, float(worst)))


def check_primes(K=100000):
    # sieve
    N = 3 * K * int(log(K + 10) + 10)
    sieve = bytearray([1]) * (N + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    ps = [p for p in range(5, N + 1, 4) if sieve[p]][:K]
    assert len(ps) == K
    s = 0.0
    for k, p in enumerate(ps, 1):
        s += log(p)
        assert p >= 4 * k + 1
        assert s >= k * log(4 * k / e) - 1e-9
    print("A.4: sum_{j<=k} log p_j >= k log(4k/e) for all k <= %d (p_j >= 4j+1 also checked): OK" % K)


if __name__ == "__main__":
    check_binary()
    check_multilevel()
    check_primes()
