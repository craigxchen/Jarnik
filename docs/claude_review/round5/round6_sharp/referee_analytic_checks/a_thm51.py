"""Referee (analytic lens): Theorem 5.1 (P3-char) step 4, re-derived constants.

(1) m_struct = prod_(q<4M) q^(a*(q)) exactly (own p'(q)), against pi(4M) log(32 M^2) and 11 M.
(2) worst-case per-prime theta-moments under the tail bounds of Lemma 3.3(3):
      small q, random lifts: Pr[(a - a*)^+ >= j] <= q^-j        ->  claim E <= 1 + 2 q^(theta-1)
      large q:               Pr[a >= a] <= min(1, 16 c / m_q) q^-(a-1) -> claim E <= 1 + 64 c q^(theta-1)
    theta = 1/(A+1) <= 1/2, levels up to A_q (X = 3 M^A); exact maximisation over tail-dominated laws.
(3) Sigma_theta = sum_(q<=X) q^(theta-1) <= X^theta/theta (exact prime sums).
(4) natural-log size of each union-bound family with the PROVED constants (b = 33, C = 1),
    as a function of log M, using tau = max(A,2) log M + 2 loglog M + log(21 b^2) (lower end):
      ternary:   M log 3 + theta log(4(M+pi/C)) - theta(tau kappa_b M/4 - 1/4) + 11 theta M + 66 Sigma_theta
      amplitude: h = 2 term: M log 5 + theta log(4(2M+pi/C)) - theta(7 tau 2M - 1/2)/2 + 11 theta M + 130 Sigma_theta + 1
      pairs:     2 log M - log 2 + log(4(2+pi/C)) + mu_M log M - (b-4) tau/4 + 1/4      (mu_M = 10 + 26/loglog M)
      off-unit:  2 log M - log 2 - theta(W/4 - 1 - 11 M) + 66 Sigma_theta,  W >= r tau
    Sigma_theta is replaced by its proved bound X^theta/theta = (A+1) 3^theta M^(A/(A+1)).
"""
import math
import numpy as np

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]

PR = primes_upto(1 << 22)
pprime = {3: 2, 5: 2, 7: 3}
k = 3
while (1 << (k + 1)) <= (1 << 22):
    qs = [int(q) for q in PR[(PR >= (1 << k)) & (PR < (1 << (k + 1)))] if q > 7]
    ts = [int(t) for t in PR[(PR >= (1 << (k - 2))) & (PR < (1 << (k - 1)))]]
    for i, q in enumerate(qs):
        pprime[q] = ts[(i * len(ts)) // len(qs)]
    k += 1
print("(1) m_struct:")
for M in (44, 1000, 10 ** 4, 10 ** 5, 10 ** 6):
    tot = 0.0
    for q in PR[(PR > 2) & (PR < 4 * M)]:
        q = int(q); a = 1
        while pprime[q] * q ** (a - 1) < M:
            a += 1
        tot += a * math.log(q)
    npi = int(((PR > 2) & (PR < 4 * M)).sum())
    print(f"   M = {M:8d}: log m_struct = {tot / M:6.3f} M;  pi(4M) log(32M^2) = {npi * math.log(32 * M * M) / M:6.3f} M  (claim <= 11 M)")

def Aq(q, X):
    a = 0
    while q ** (a + 1) <= X:
        a += 1
    return a

w_small = w_large = 0.0
for A in (1, 2, 3):
    theta = 1.0 / (A + 1)
    for M in (44, 10 ** 3, 10 ** 5):
        X = 3 * M ** A
        for q in (3, 5, 7, 11, 101, 4 * M + 1, 10 * M + 1, int(X ** 0.5) + 1, X // 3):
            if q < 3 or q > X:
                continue
            Am = Aq(q, X)
            # small-q random lifts above a* (take a* = 1 as worst case: most random levels)
            E = 1.0
            for j in range(1, Am):
                E += q ** (-j) * (q ** (theta * j) - q ** (theta * (j - 1)))
            w_small = max(w_small, (E - 1) / (2 * q ** (theta - 1)))
            if q > 4 * M:
                for c in (1, 2, 5, 50):
                    m = q - 1 if q % 4 == 1 else q + 1
                    E = 1.0
                    for a in range(1, Am + 1):
                        t = min(1.0, 16 * c / m) * q ** (-(a - 1))
                        E += t * (q ** (theta * a) - q ** (theta * (a - 1)))
                    w_large = max(w_large, (E - 1) / (64 * c * q ** (theta - 1)))
print(f"(2) worst (E-1)/(2 q^(theta-1)) small-q lifts: {w_small:.3f};  worst (E-1)/(64 c q^(theta-1)) large q: {w_large:.3f}  (<= 1 required)")

worst = 0.0
for A in (1, 2, 3):
    theta = 1.0 / (A + 1)
    for M in (44, 300, 3000):
        X = 3 * M ** A
        if X > 3 * 10 ** 7:
            continue
        P = primes_upto(X)
        P = P[P > 2].astype(float)
        S = (P ** (theta - 1)).sum()
        worst = max(worst, S / (X ** theta / theta))
print(f"(3) max Sigma_theta / (X^theta/theta) = {worst:.3f}  (<= 1 required)")

print("(4) ln(family total) with proved constants, b = 33, C = 1:")
b, C = 33, 1.0
kap = min(1.0, b / 16 - 2)
for A in (1, 2):
    theta = 1.0 / (A + 1)
    row = []
    for L in (50, 100, 200, 400, 450, 1000, 1e4, 3.5e4, 1e5):
        M = math.exp(min(L, 700))   # only used where it matters numerically via logs below
        ll = math.log(L)
        tau = max(A, 2) * L + 2 * ll + math.log(21 * b * b)
        logSig = math.log(A + 1) + theta * math.log(3) + (A / (A + 1)) * L
        def big(x):   # e^x with overflow guard (x in log space)
            return x
        # ternary: exponent = M * [log 3 - theta(tau kap/4 - 11)] + small; report sign of bracket and the Sigma term relative
        brk_t = math.log(3) - theta * (tau * kap / 4 - 11)
        brk_a = math.log(5) / 2 - theta * (7 * tau - 11 / 2)    # per (h M) for h = 2, conservative
        mu = 10 + 26 / ll
        pairs = 2 * L - math.log(2) + math.log(4 * (2 + math.pi / C)) + mu * L - (b - 4) * tau / 4 + 0.25
        # Sigma_theta/M = e^(logSig - L)
        sig_over_M = math.exp(logSig - L)
        row.append((L, brk_t, 66 * sig_over_M, brk_a, pairs))
    for (L, bt, sg, ba, pa) in row:
        print(f"   A = {A}, log M = {L:7.0f}: ternary exponent/M = {bt + sg:9.2f}  (log3 - theta(tau kappa_b/4 - 11) + 66 Sigma/M);"
              f"  amplitude-2 exponent/(2M) = {ba:9.2f};  pairs ln = {pa:10.1f}")
