"""Necessary condition of Prop 5.1 for cheap spans at m' = ceil(4w) (nearstrata.md Prop 5.3).
For every cheap family F (Theorem 3.4): s = sum ell <= b = floor((m'-1)/(2j)), |E| = e <= e_max(s) = floor((m'-1)/2) - j s,
q_F <= 2^s, P(L_F) >= n'(n'-q)/(2q), n' = n - e.   R(s,e) := 2^{s+1} N_2 / ((n-e)(n-e-2^s)) >= N_2/P(L_F).
 Condition A (any true transversal order o <= m'):   max_{s,e} m' R(s,e) <= d.
 Condition B (proved order lower bound o >= m' - d_F >= m' - 2js - 2e):  max_{s,e} (m' - 2js - 2e) R(s,e) <= d.
R(s,e) is increasing in e, so for A the max over e is at e = e_max(s); for B all e are scanned (numpy).
Scan: every n in [100, 20000], then geometric step 1.0005 up to 4*10^6; at each reported threshold the next 200 n are
re-checked, and condition A is re-evaluated in exact rational arithmetic at threshold + 1.
"""
import math
import numpy as np
from fractions import Fraction as Fr


def params(n, r):
    d = n - 1
    w = math.sqrt(8 * d * math.log(4 * r))
    return d, w, math.floor(w / 2), math.ceil(4 * w)


def conds_float(n, r):
    d, w, j, mp = params(n, r)
    if w < 10 or 2 * j > n:
        return None
    N2 = n * (n - 1) / 2
    b = (mp - 1) // (2 * j)
    A = B = 0.0
    for s in range(b + 1):
        emax = (mp - 1) // 2 - j * s
        if emax < 0:
            continue
        q = 2.0 ** s
        e = np.arange(emax + 1, dtype=float)
        R = 2 * q * N2 / ((n - e) * (n - e - q))
        A = max(A, mp * R[-1])
        B = max(B, float(np.max(np.maximum(mp - 2 * j * s - 2 * e, 0) * R)))
    return A / d, B / d


def condA_exact(n, r):
    d, w, j, mp = params(n, r)
    N2 = n * (n - 1) // 2
    b = (mp - 1) // (2 * j)
    A = Fr(0)
    for s in range(b + 1):
        emax = (mp - 1) // 2 - j * s
        if emax < 0:
            continue
        q = 2 ** s
        A = max(A, mp * Fr(2 * q * N2, (n - emax) * (n - emax - q)))
    return A / d


def scan(rf, hi):
    lastA = lastB = None
    n = 100
    while n <= hi:
        c = conds_float(n, rf(n))
        if c is None:
            lastA = lastB = n
        else:
            if c[0] > 1:
                lastA = n
            if c[1] > 1:
                lastB = n
        n = n + 1 if n < 20000 else int(n * 1.0005) + 1
    return lastA, lastB


if __name__ == '__main__':
    for name, rf in [('r = n', lambda n: n), ('r = n ln n', lambda n: n * math.log(n)), ('r = n^2', lambda n: n * n)]:
        a, b = scan(rf, 4 * 10**6)
        okA = all(conds_float(n, rf(n))[0] <= 1 for n in range(a + 1, a + 201))
        okB = all(conds_float(n, rf(n))[1] <= 1 for n in range(b + 1, b + 201))
        exA = float(condA_exact(a + 1, rf(a + 1)))
        print(f"{name}: condition A last failing n = {a} (next 200 pass: {okA}; exact A-ratio at {a+1}: {exA:.6f}); "
              f"condition B last failing n = {b} (next 200 pass: {okB})")
    print()
    print("n (r = n): b, A-ratio = max m' N2/(P d), B-ratio, m', 8w, 64w")
    for n in [3000, 5185, 10**4, 10**5, 4 * 10**5, 5 * 10**5, 10**6, 2 * 10**6]:
        d, w, j, mp = params(n, n)
        A, B = conds_float(n, n)
        print(n, (mp - 1) // (2 * j), round(A, 4), round(B, 4), mp, round(8 * w, 1), round(64 * w, 1))
