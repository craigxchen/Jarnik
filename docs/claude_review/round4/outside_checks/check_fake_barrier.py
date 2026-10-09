"""Exact checks for the combinatorial core of Theorem A.3 (fake-circle barrier) in outside.md.

Profile: Hadamard H of order M (first column all ones), labels a = 1..M-1, b copies per label;
copies 0,1 of every label are flipped (each at one row), copies 2,3 are their siblings.
Every row receives 1 or 2 flips and one designated flip (twin = designated flip + its sibling).
S_{x,j} = H(x, a(j)) * (-1)^[j flipped at row x].   For zero-sum c in Z^M (n = ||c||_1):
   T_a = sum_x c_x H(x,a),  F(c) = sum_{a>=1} |T_a|,  G(c) = sum_j |c^T S_j| - b(M-1).
Equal-weight forms of the inequalities proved in the text (weights = 1):
   (F)   F(c) >= M
   (B1u) G >= b (F - M + 1) - 4 n
   (B2u) G >= 2 n + (b - 4)(F - M + 1) - 4 (M - 1)
   (fin) with b = 6 M:   G/2 >= n + M      (equal-weight form of (B2) at b = 6M)
Also: S has rank M (saturation: every rational character with integral image is integral),
and the pair characters c = e_x - e_y have G >= b - 8 (two flips per row at worst).
All arithmetic is exact (Python integers).
"""
import random
from itertools import combinations, product
import numpy as np

random.seed(5)


def sylvester(m):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < m:
        H = np.block([[H, H], [H, -H]])
    return H


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley1(q):
    # q = 3 mod 4 prime; order q+1
    Q = np.array([[legendre(j - i, q) for j in range(q)] for i in range(q)], dtype=np.int64)
    n = q + 1
    S = np.zeros((n, n), dtype=np.int64)
    S[0, 1:] = 1
    S[1:, 0] = -1
    S[1:, 1:] = Q
    H = np.eye(n, dtype=np.int64) + S
    return H


def normalize(H):
    H = H * H[:, [0]]  # row signs so first column is all ones
    M = H.shape[0]
    assert (H.T @ H == M * np.eye(M, dtype=np.int64)).all()
    assert (H[:, 0] == 1).all()
    return H


def build(H, b):
    M = H.shape[0]
    cols = []  # (label, copy)
    for a in range(1, M):
        for t in range(b):
            cols.append((a, t))
    r = len(cols)
    S = np.zeros((M, r), dtype=np.int64)
    flip_row = {}
    # flipped columns: copies 0 and 1 of each label; assign rows cyclically
    flipped = [(a, t) for t in (0, 1) for a in range(1, M)]
    designated = {}
    for i, col in enumerate(flipped):
        x = i % M
        flip_row[col] = x
        if i < M:
            designated[x] = col
    idx = {col: j for j, col in enumerate(cols)}
    for j, (a, t) in enumerate(cols):
        S[:, j] = H[:, a]
        if (a, t) in flip_row:
            S[flip_row[(a, t)], j] *= -1
    flips_per_row = [sum(1 for col, x in flip_row.items() if x == y) for y in range(M)]
    assert min(flips_per_row) >= 1 and max(flips_per_row) <= 2
    return S, flip_row, designated, idx


def quantities(H, S, b, c):
    M = H.shape[0]
    T = c @ H
    assert T[0] == 0
    F = int(np.abs(T[1:]).sum())
    G = int(np.abs(c @ S).sum()) - b * (M - 1)
    return F, G


def test(H, b, name, nrand=3000):
    M = H.shape[0]
    S, flip_row, designated, idx = build(H, b)
    assert np.linalg.matrix_rank(S.astype(float)) == M
    # twin structure: designated flipped column minus its sibling = -2 H(x,a) e_x
    for x, (a, t) in designated.items():
        d = S[:, idx[(a, t)]] - S[:, idx[(a, t + 2)]]
        e = np.zeros(M, dtype=np.int64)
        e[x] = -2 * H[x, a]
        assert (d == e).all()
    cs = []
    # all zero-sum c with ||c||_1 <= 4
    for k in (2, 3, 4):
        for supp in combinations(range(M), min(k, M)):
            for vals in product(*[range(-4, 5)] * len(supp)):
                if 0 in vals or sum(vals) != 0 or sum(abs(v) for v in vals) > 4:
                    continue
                c = np.zeros(M, dtype=np.int64)
                c[list(supp)] = vals
                cs.append(c)
        if len(cs) > 60000:
            break
    # columns, two-column half sums, random vectors of all sizes
    for a in range(1, M):
        cs.append(H[:, a].copy())
        for a2 in range(a + 1, M):
            for s in (1, -1):
                cs.append((H[:, a] + s * H[:, a2]) // 2)
    for _ in range(nrand):
        n = random.randint(2, 4 * M)
        c = np.zeros(M, dtype=np.int64)
        for _ in range(n // 2):
            i, j = random.sample(range(M), 2)
            c[i] += 1
            c[j] -= 1
        if not c.any():
            continue
        cs.append(c)
    worst_fin = None
    worst_pair = None
    for c in cs:
        n = int(np.abs(c).sum())
        F, G = quantities(H, S, b, c)
        assert F >= M, (name, c)
        assert G >= b * (F - M + 1) - 4 * n, (name, c)
        assert G >= 2 * n + (b - 4) * (F - M + 1) - 4 * (M - 1), (name, c)
        if b == 6 * M:
            slack = G / 2 - (n + M)
            assert slack >= 0, (name, c, G)
            worst_fin = slack if worst_fin is None else min(worst_fin, slack)
        if n == 2:
            assert G >= b - 8
            worst_pair = G if worst_pair is None else min(worst_pair, G)
    print("%-14s M=%3d b=%4d: %6d characters checked; (F),(B1u),(B2u) OK; min pair G=%d (b-8=%d)%s"
          % (name, M, b, len(cs), worst_pair, b - 8,
             "" if worst_fin is None else "; (fin) min slack %.1f" % worst_fin))


def main():
    for m in (8, 16, 32):
        H = normalize(sylvester(m))
        test(H, 6, "Sylvester", 1000)
        test(H, 6 * m, "Sylvester", 1000)
    for q in (11, 19, 23):
        H = normalize(paley1(q))
        test(H, 6, "Paley-I", 1000)
        test(H, 6 * (q + 1), "Paley-I", 1000)


if __name__ == "__main__":
    main()
