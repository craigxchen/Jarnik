#!/usr/bin/env python3
"""Referee check of Prop. D.1 (p3lower.md): clique-partition sign matrices, exact integers.
Also exhibits that the lower bound r >= lambda_max(k) needs k zero-diagonal ('[same group]' must exclude x=y)."""
import numpy as np, itertools, math
from sympy import Matrix

def sylv(n):
    H = np.array([[1]])
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    assert H.shape[0] == n
    return H

fails = 0
cases = [(2, 2, 1), (4, 2, 5), (4, 4, 5), (8, 2, 9), (8, 4, 1), (2, 8, 3), (4, 8, 13), (16, 4, 17), (8, 16, 9), (32, 4, 97)]
for g, t, L in cases:
    m = g - 1 + L
    if m & (m - 1):  # Sylvester orders only
        print("skip", g, t, L); continue
    Hg, Ct, Gm = sylv(g), sylv(t), sylv(m)
    M = g * t
    rows = []
    for a in range(g):
        for i in range(t):
            rows.append(np.concatenate([Hg[a, 1:], np.kron(Ct[i, 1:], Gm[a, :m])]))
    S = np.array(rows)
    r = S.shape[1]
    G = S @ S.T
    # expected Gram
    ok = r == M - t + (t - 1) * L
    for x, (a, i) in enumerate(itertools.product(range(g), range(t))):
        for y, (b, j) in enumerate(itertools.product(range(g), range(t))):
            exp = r if x == y else (-L if a == b else -1)
            ok &= G[x, y] == exp
    ok &= all(len(set(S[:, c])) == 2 for c in range(r))  # nonconstant columns
    k = np.ones((M, M), dtype=int) - np.eye(M, dtype=int)
    for x, (a, i) in enumerate(itertools.product(range(g), range(t))):
        for y, (b, j) in enumerate(itertools.product(range(g), range(t))):
            if x != y and a == b: k[x, y] = L
    lam = max(np.linalg.eigvalsh(k.astype(float)))
    ok &= abs(lam - r) < 1e-6
    fails += (not ok)
    print(f"g={g:3d} t={t:3d} L={L:3d} M={M:4d}: r={r:5d}, lambda_max(k)={lam:9.3f}, uniform (M-1)L={(M-1)*L:6d}  {'OK' if ok else 'FAIL'}")
print("ALL D.1 CHECKS PASSED" if fails == 0 else f"FAILURES {fails}")
