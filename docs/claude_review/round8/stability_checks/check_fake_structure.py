#!/usr/bin/env python3
"""The local structure theorems on a refereed fake circle (round5/construct.md Theorem 2, exact instance
M = 12, b = 8, C = 1, primes = the first r = 100 primes = 1 mod 4 above 10^7).

The profile is rebuilt with the code of round5/checks/lib.py (copied here: paley1, normalize,
flipped_profile).  For the angle model only the profile and the prime weights matter for the local
statements; Proposition 7.1 of stability.md says every conclusion of Sections 2-5 must hold.  We print:
  * W, A = W/(M log M), max log p against the smoothness bound 2W/M - 3 log(M/2) + K0 and b(s);
  * the sub-cluster gcd bound gamma_B <= W/s - 8F*(s)/s^2 (C = 1) for all 4083 subsets B, |B| >= 2;
  * omega = r, the number of distinct cut partitions, cut slack Q, the margin range, the minimum
    Hamming distance between rows, and the minimum character support over pairs and Hadamard columns.
Weights are the exact prime logarithms (mpmath, 60 digits).
"""
from itertools import combinations
import numpy as np
import mpmath
from stab_lib import Fstar, K0

mpmath.mp.dps = 60


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley1(q):
    Q = np.array([[legendre(j - i, q) for j in range(q)] for i in range(q)], dtype=np.int64)
    n = q + 1
    S = np.zeros((n, n), dtype=np.int64)
    S[0, 1:] = 1
    S[1:, 0] = -1
    S[1:, 1:] = Q
    return np.eye(n, dtype=np.int64) + S


def normalize(H):
    H = H * H[:, [0]]
    M = H.shape[0]
    assert (H.T @ H == M * np.eye(M, dtype=np.int64)).all()
    return H


def flipped_profile(H, b, assign):
    M = H.shape[0]
    cols = []
    for a in range(1, M):
        for _ in range(b):
            cols.append(H[:, a].copy())
    for x in range(M):
        col = H[:, assign[x]].copy()
        col[x] = -col[x]
        cols.append(col)
    return np.array(cols, dtype=np.int64).T


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def main():
    q, b = 11, 8
    M = q + 1
    H = normalize(paley1(q))
    assign = [1] + list(range(1, q + 1))
    S = flipped_profile(H, b, assign)
    r = S.shape[1]
    ps = []
    n = 10 ** 7 + 1
    while len(ps) < r:
        if n % 4 == 1 and is_prime(n):
            ps.append(n)
        n += 1
    w = [mpmath.log(p) for p in ps]
    W = sum(w)
    print(f"fake instance: Paley q={q}, M={M}, b={b}, r={r}, primes {ps[0]}..{ps[-1]}")
    print(f"  W = {float(W):.2f},  A = W/(M log M) = {float(W / (M * mpmath.log(M))):.2f}")
    bnd = 2 * W / M - 3 * mpmath.log(M / 2) + K0
    print(f"  max log p = {float(max(w)):.3f}  vs smoothness bound 2W/M - 3log(M/2) + K0 = {float(bnd):.2f}")
    # sub-cluster gcd lemma on all subsets
    worst = None
    nfail = 0
    cnt = 0
    for s in range(2, M + 1):
        bs = W / s - 8 * Fstar(s) / (s * s)
        for B in combinations(range(M), s):
            sub = S[list(B), :]
            const = np.all(sub == sub[0], axis=0)
            gam = sum(w[j] for j in range(r) if const[j])
            cnt += 1
            sl = gam - bs
            if sl > 0:
                nfail += 1
            if worst is None or sl > worst[0]:
                worst = (float(sl), s, int(const.sum()))
    print(f"  gcd lemma on {cnt} subsets: failures {nfail}; max gamma_B - b(s) = {worst[0]:.2f} "
          f"(s = {worst[1]}, {worst[2]} constant cuts)")
    # partitions, cut slack, margins, Hamming distance
    parts = set()
    for j in range(r):
        col = tuple(S[:, j] if S[0, j] == 1 else -S[:, j])
        parts.add(col)
    Q = sum(w[j] * mpmath.mpf(int(S[:, j].sum())) ** 2 / 4 for j in range(r))
    print(f"  omega = r = {r} (M/4 = {M / 4}),  distinct cut partitions = {len(parts)} (> M/2 = {M // 2})")
    print(f"  cut slack Q = {float(Q):.2f}  vs  MW/4 = {float(M * W / 4):.1f}")
    margins = []
    hmin = None
    for x in range(M):
        for y in range(x):
            d = sum(w[j] for j in range(r) if S[x, j] != S[y, j])
            margins.append(d - W / 2)
            h = int((S[x] != S[y]).sum())
            hmin = h if hmin is None else min(hmin, h)
    K = M * (M - 1) // 2
    print(f"  margins m_xy: min {float(min(margins)):.2f}, mean {float(sum(margins) / K):.2f}, "
          f"max {float(max(margins)):.2f};  W/(2(M-1)) = {float(W / (2 * (M - 1))):.2f};  log M = {float(mpmath.log(M)):.2f}")
    print(f"  min Hamming distance between rows = {hmin}  (> M/4 = {M / 4})")
    # character supports: pairs and Hadamard columns (c = h_a)
    sup_pairs = min(int(((S[x] - S[y]) != 0).sum()) for x in range(M) for y in range(x))
    sup_cols = min(int(((S.T @ H[:, a]) != 0).sum()) for a in range(1, M))
    print(f"  min support of pair characters = {sup_pairs};  of Hadamard-column characters = {sup_cols}")
    print("ALL CHECKS PASSED:", nfail == 0 and len(parts) > M / 2 and hmin > M / 4)


if __name__ == "__main__":
    main()
