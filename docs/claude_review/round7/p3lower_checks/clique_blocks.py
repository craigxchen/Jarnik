"""Exact check of Proposition D.1 (clique-partition blocks).

Points are split into g groups of size t (M = g t).  We want +-1 rows s_x (equal weights) with
    <s_x, s_y> = -L  (same group),   <s_x, s_y> = -1  (different groups).
Construction (Sylvester Hadamard matrices, all orders powers of 2):
    H_g : Hadamard of order g  (row per group; nonconstant columns give -1 between groups)
    C_t : Hadamard of order t  (row per member; nonconstant columns give -1 inside a group)
    G_m : Hadamard of order m >= g, m = g - 1 + L  (orthogonal rows, one per group)
    s_(group a, member i) = ( H_g[a, 1:],  C_t[i, 1:] (x) G_m[a, :] )
Then r = (g - 1) + (t - 1) m = M - t + (t - 1) L = lam_max(k) for k = (J - I) + (L - 1)[same group],
i.e. the Perron lower bound r >= lam_max(k) is attained.  t = 2 is the matching case.
Integer arithmetic throughout; this is an exact verification of the construction for the listed
parameters (the proof in the write-up is general).
"""
import numpy as np

def sylvester(n):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    assert H.shape[0] == n
    return H

def build(g, t, L):
    m = g - 1 + L
    Hg, Ct, Gm = sylvester(g), sylvester(t), sylvester(m)
    rows = []
    for a in range(g):
        for i in range(t):
            rows.append(np.concatenate([Hg[a, 1:], np.kron(Ct[i, 1:], Gm[a, :])]))
    return np.array(rows)

def check(g, t, L):
    S = build(g, t, L)
    M, r = S.shape
    G = S @ S.T
    ok = True
    for x in range(M):
        for y in range(M):
            if x == y:
                ok &= G[x, y] == r
            elif x // t == y // t:
                ok &= G[x, y] == -L
            else:
                ok &= G[x, y] == -1
    # nonconstant columns
    ok &= bool(np.all(np.abs(S.sum(axis=0)) < M))
    k = np.ones((M, M), dtype=float) - np.eye(M)
    for x in range(M):
        for y in range(M):
            if x != y and x // t == y // t:
                k[x, y] = L
    lam = max(np.linalg.eigvalsh(k))
    uniform_r = (M - 1) * L  # cost if every pair gets correlation -L (Hadamard copies)
    print(f"g={g:3d} t={t:2d} L={L:3d}  M={M:4d} r={r:5d} lam_max(k)={lam:9.3f} "
          f"(M-t+(t-1)L={M - t + (t - 1) * L})  uniform cost (M-1)L={uniform_r:6d}  exact Gram: {ok}")
    return ok and abs(lam - r) < 1e-6

allok = True
for (g, t, L) in [(2, 2, 3), (4, 2, 5), (8, 2, 9), (16, 2, 17), (16, 2, 49), (32, 2, 33),
                  (4, 4, 5), (8, 4, 9), (8, 8, 25), (16, 4, 17), (4, 16, 13), (32, 4, 97)]:
    allok &= check(g, t, L)
print("ALL CLIQUE-BLOCK CHECKS PASSED" if allok else "SOME CHECK FAILED")
