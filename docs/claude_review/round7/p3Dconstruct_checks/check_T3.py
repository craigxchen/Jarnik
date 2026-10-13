"""Theorem 3 checks: Hadamard uncertainty (step 3), row-image identity (step 1), and two explicit
G-coherent designs on the sharp.md section-4 profile showing the dichotomy (small columns => large pairs;
small pairs => degenerate label sums).  Exact integer/rational arithmetic (fractions) where it matters."""
import itertools, random
from fractions import Fraction
import numpy as np
from gauss import legendre, primes_upto

def paley_normalised(q0):
    M = q0 + 1
    S = np.zeros((M, M), dtype=int)       # index 0 = infinity, 1+a = a in F_q0
    for a in range(q0):
        S[0, 1+a] = 1; S[1+a, 0] = -1
        for b in range(q0):
            S[1+a, 1+b] = legendre(b - a, q0)
    H = np.eye(M, dtype=int) + S
    H = H * H[:, [0]]                     # make column 0 (infinity) all ones
    assert (H @ H.T == M*np.eye(M, dtype=int)).all() and (H[:, 0] == 1).all()
    return H

def uncertainty(q0, rng):
    H = paley_normalised(q0); M = q0+1; Hp = H[:, 1:]
    worst = None
    cands = []
    for a in range(M-1):
        for c in (1, 2, -3):
            psi = np.zeros(M-1, dtype=int); psi[a] = c; cands.append(psi)
    for _ in range(3000):
        psi = np.zeros(M-1, dtype=int)
        k = rng.choice([2, 3, 5, 10, M-1])
        for a in rng.sample(range(M-1), k): psi[a] = rng.choice([-2, -1, 1, 2])
        cands.append(psi)
    for psi in cands:
        y = Hp @ psi / 2.0
        val = np.abs(y - np.median(y)).sum()          # median minimises the L1 deviation
        if worst is None or val < worst: worst = val
    return M, worst

def profile(q0, b, astar=0):
    H = paley_normalised(q0); M = q0+1; col = (1 + H[:, 1:]) // 2   # 0/1 label columns, labels 0..q0-1
    cols, idx = [], {}
    for a in range(q0):
        for i in range(b):
            idx[(a, i)] = len(cols); cols.append(col[:, a].copy())
    A = np.array(cols).T.copy()            # M x r
    twins = {}
    for x in range(M):
        if x == 0: f, s_ = (astar, 2), (astar, 3)
        else: f, s_ = (x-1, 0), (x-1, 1)
        j = idx[f]; eps = 1 if A[x, j] == 0 else -1
        A[x, j] += eps; twins[x] = (idx[f], idx[s_], eps)
    return A, idx, twins

def in_L(A, v):
    """v in L = {A^T c : sum c = 0} over Q (L is saturated by sharp Lemma 4.1)."""
    M = A.shape[0]
    Mat = np.vstack([A.T, np.ones((1, M))])
    rhs = np.concatenate([v, [0]])
    c, res, rk, _ = np.linalg.lstsq(Mat.astype(float), rhs.astype(float), rcond=None)
    return np.abs(Mat @ c - rhs).max() < 1e-8

def pair_stats(G, A):
    M = A.shape[0]; P = G @ A.T          # D x M ; column x = G a_x
    Hs = [np.abs(P[:, x] - P[:, y]).sum() for x in range(M) for y in range(x+1, M)]
    return np.mean(Hs), max(Hs)

def main():
    rng = random.Random(7)
    for q0 in (43, 47, 59):
        M, w = uncertainty(q0, rng)
        print(f"uncertainty q0={q0}: min ||H'psi/2 - median||_1 = {w:.1f}  vs M/4 = {M/4}")
        assert w >= M/4 - 1e-9
    q0, b = 43, 6
    A, idx, twins = profile(q0, b); M, r = A.shape
    # step 1 identity on random G rows
    for _ in range(50):
        g = np.array([rng.randrange(-2, 3) for _ in range(r)])
        psi = np.array([sum(g[idx[(a, i)]] for i in range(b)) for a in range(q0)])
        H = paley_normalised(q0); Hp = H[:, 1:]
        phi = np.zeros(M)
        for x, (jf, js, eps) in twins.items(): phi[x] = eps * g[jf]
        F = Hp @ psi / 2 + phi
        d = A @ g - F
        assert np.allclose(d, d[0])
    print("row-image identity (step 1): ok on 50 random rows")
    # design 1: G = identity (each column its own Gaussian prime): C = 1
    G1 = np.eye(r, dtype=int)
    m1, x1 = pair_stats(G1, A)
    print(f"design G=I: C=1, mean pair height {m1:.1f}, max {x1}  (r = {r}, M = {M})")
    # design 2: rows = within-label path differences (+ twin rows automatically among them)
    rows = []
    for a in range(q0):
        for i in range(b-1):
            g = np.zeros(r, dtype=int); g[idx[(a, i)]] = 1; g[idx[(a, i+1)]] = -1; rows.append(g)
    G2 = np.array(rows)
    C2 = np.abs(G2).sum(axis=0).max()
    m2, x2 = pair_stats(G2, A)
    # degenerate monomial: a label sum
    v = np.zeros(r, dtype=int)
    for i in range(b): v[idx[(5, i)]] = 1
    print(f"design G=paths: C={C2}, mean pair height {m2:.2f}, max {x2}; label sum s_5: Gv=0 {not (G2@v).any()},"
          f" in L {in_L(A, v)}, weight {v.sum()}")
    assert not (G2 @ v).any() and not in_L(A, v)
    # pigeonhole demo with design 2 plus a few dense rows: count how many rows see label classes
    print("ALL T3 CHECKS PASSED")

if __name__ == "__main__":
    main()
