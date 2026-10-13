#!/usr/bin/env python3
"""Referee check of the E(n,Q) facts and of the 'factor 2' claims around Prop. B.3 (p3lower.md).
1. exhaustive: E formula, n^2/(2Q)-n/2 <= E <= n^2/(2Q), superadditivity, B.3(3) chain.
2. per-modulus ratio max_{n0} [E(n0,m/2)+E(M-n0,m/2)] / E(M,m): NOT bounded by 2 (Sec. 0 table claim).
3. aggregated B.3(4) main terms (8/M) sum_Q log q * count, parity-free vs full bias, X = 3M:
   compare with the claimed (4-o(1)) M log M and (8-o(1)) M log M (Sec. 6 table)."""
import math, itertools
from sympy import primerange

def E(n, Q):
    b, r = divmod(n, Q)
    return Q * b * (b - 1) // 2 + r * b

def Ebrute(n, Q):
    # least collisions: balanced distribution (convexity); brute force over compositions for tiny n
    best = None
    def rec(k, left, acc):
        nonlocal best
        if k == Q - 1:
            v = acc + left * (left - 1) // 2
            best = v if best is None or v < best else best
            return
        for c in range(left + 1):
            rec(k + 1, left - c, acc + c * (c - 1) // 2)
    rec(0, n, 0)
    return best

fails = 0
for Q in range(1, 6):
    for n in range(0, 11):
        if E(n, Q) != Ebrute(n, Q): fails += 1
for Q in range(1, 13):
    for n in range(0, 60):
        if not (n * n / (2 * Q) - n / 2 - 1e-12 <= E(n, Q) <= n * n / (2 * Q) + 1e-12): fails += 1
    for n0 in range(0, 40):
        for n1 in range(0, 40):
            if E(n0, Q) + E(n1, Q) > E(n0 + n1, Q): fails += 1
for m in range(4, 25, 4):
    for n0 in range(0, 40):
        for n1 in range(0, 40):
            Mm = n0 + n1
            mid = E(n0, m // 2) + E(n1, m // 2)
            if not (E(Mm, m) <= mid <= E(Mm, m // 2) <= Mm * Mm / m + 1e-12): fails += 1
            beta = (n0 - n1) / Mm if Mm else 0
            if mid < Mm * Mm / (2 * m) * (1 + beta ** 2) - Mm / 2 - 1e-9: fails += 1
print("E facts:", "OK" if fails == 0 else f"FAIL {fails}")

print("per-modulus parity gain ratios (worst n0):")
for (M, m) in [(10, 8), (10, 12), (10, 16), (40, 36), (40, 60), (100, 12), (100, 120)]:
    pf = E(M, m)
    best = max(E(n0, m // 2) + E(M - n0, m // 2) for n0 in range(M + 1))
    print(f"  M={M:4d} m={m:4d}: parity-free E(M,m)={pf:5d}, worst parity count={best:5d}, ratio={'inf' if pf == 0 else round(best/pf, 3)}")

print("aggregated B.3(4) main terms / (M log M), X = 3M, all odd q <= X assumed coprime to N:")
for M in [100, 1000, 10000, 100000]:
    X = 3 * M
    s_free = 0.0; s_bias = 0.0; s_bal = 0.0
    for q in primerange(3, X + 1):
        a = 1
        while q ** a <= X:
            m = (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)
            lq = math.log(q)
            s_free += lq * E(M, m)
            s_bias += lq * E(M, m // 2)
            s_bal += lq * (E(M // 2, m // 2) + E(M - M // 2, m // 2))
            a += 1
    f = 8 / M
    L = M * math.log(M)
    print(f"  M={M:6d}: parity-free {f*s_free/L:.3f}, balanced parity {f*s_bal/L:.3f}, full bias {f*s_bias/L:.3f}; aggregated ratio bias/free = {s_bias/s_free:.3f}")
print("ALL E CHECKS DONE" if fails == 0 else f"FAILURES {fails}")
