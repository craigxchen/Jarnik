#!/usr/bin/env python3
"""Referee check (block fakedefs): where Astra's VIII.3 profile sits relative to the record's window and to
round-21 Theorem A, Corollary A.2 and judge J1.

Part 1 (actual primes, moderate q): P = q^d; find the shortest multiplicative interval [p_a, p_b] in [P, 2P]
containing q-1 split primes; sigma = log(p_b/p_a). Compare with the J1 equal-weight width 1/(4r), r = q-1, and
with the J1 threshold 16 r (r-1) e^{-1/(4r)}. Random permutation of these primes on the q-1 columns of the Paley
profile at q = 383 (S = class A, rho = 1 limit), exact pair margins w(K_i - K_i') - W/2 against 2 log p_a.
Part 2 (asymptotic bookkeeping, no primes): exponents log P_+/log m, W/(m log m), Lambda, Theorem A exponent.
"""
import math, random
import numpy as np
from sympy import primerange


def chi(x, q):
    x %= q
    if x == 0:
        return 0
    return 1 if pow(x, (q - 1) // 2, q) == 1 else -1


def shortest_window(q, d):
    P = int(round(q ** d))
    ps = [p for p in primerange(P, 2 * P + 1) if p % 4 == 1]
    r = q - 1
    best = None
    for i in range(len(ps) - r + 1):
        s = math.log(ps[i + r - 1] / ps[i])
        if best is None or s < best[0]:
            best = (s, i)
    s, i = best
    return P, len(ps), ps[i:i + r], s


print("Part 1: shortest window of q-1 split primes inside [P, 2P], P = q^d")
for q, d in [(103, 1.99), (383, 1.99), (1019, 1.99), (383, 23 / 12), (1019, 23 / 12)]:
    P, npr, win, sigma = shortest_window(q, d)
    r = q - 1
    j1 = 16 * r * (r - 1) * math.exp(-1 / (4 * r))
    print(f" q={q} d={d:.4f} P={P} split primes in [P,2P]={npr} sigma={sigma:.3e} 1/(4r)={1/(4*r):.3e} "
          f"sigma*4r={sigma*4*r:.1f}  J1 threshold 16r(r-1)e^(-1/4r)={j1:.3e} vs p_min={win[0]} "
          f"(J1 window possible: {win[0] >= j1}); sigma*sqrt(q log q)={sigma*math.sqrt(q*math.log(q)):.3f}")

# pair margins at q = 383 with random permutation of the window primes
q, d = 383, 1.99
P, npr, win, sigma = shortest_window(q, d)
K = np.array([[chi(j - i, q) for j in range(q)] for i in range(q)], dtype=np.int64)
A = [i for i in range(q) if K[i, 0] == 1]
KS = K[np.ix_(A, list(range(1, q)))]
logs = np.array([math.log(p) for p in win])
rng = random.Random(7)
mins = []
for trial in range(5):
    perm = logs.copy()
    rng.shuffle(perm)
    W = 2 * perm.sum()
    worst = 1e9
    for x in range(len(A)):
        dif = np.abs(KS[x + 1:] - KS[x])          # rows y > x
        marg = dif @ perm - W / 2
        if len(marg):
            worst = min(worst, marg.min())
    mins.append(worst)
print(f"Part 1b: q=383, |S|=|A|={len(A)}, sigma={sigma:.3e}: min pair margin over 5 random permutations "
      f"{[round(t,3) for t in mins]} vs 2 log p_min = {2*math.log(win[0]):.3f}; "
      f"equal-weight worst case 2 log p_min - (q-1) sigma = {2*math.log(win[0]) - (q-1)*sigma:.3f}")

print("Part 2: asymptotic bookkeeping (m = rho q/2, r = q-1, P_+ <= 2 q^d, W = 2 (q-1) log P(1+o(1)))")
for d, rho in [(1.99, 0.995), (23 / 12, 0.999)]:
    for k in [3, 6, 12, 30, 100, 200]:
        q = 10.0 ** k
        m = rho * q / 2
        logP = d * math.log(q)
        W = 2 * (q - 1) * logP
        lpm = (math.log(2) + logP) / math.log(m)
        coef = W / (m * math.log(m))
        m_free = (q - 1) - (m - 1)
        thmA = 2.1735 * m_free * math.exp(-0.5 * (math.log(2) + logP))
        Lam_low = ((q + 1) / 2) * math.exp(-0.5 * (math.log(2) + logP))
        print(f" d={d:.4f} rho={rho} q=1e{k}: log P_+/log m={lpm:.4f} W/(m log m)={coef:.3f} "
              f"(limit 4d/rho={4*d/rho:.3f}) Lambda>={Lam_low:.3g} ThmA: Haar prop <= exp(-{thmA:.3g}) "
              f"m_free/m={m_free/m:.3f}")
