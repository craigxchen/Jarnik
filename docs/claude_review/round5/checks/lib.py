"""Shared helpers for round 5 (construction of cheap angle models).

Conventions (see construct.md):
  * H : normalised Hadamard matrix of order M (first column all ones), labels a = 1..M-1.
  * A profile is a +-1 matrix S (M x r); row x of the exponent profile is a_x = (1 + S_x)/2.
  * For an integral zero-sum c, y = S^T c (all entries even), and with equal weights tau
    V_c - W/2 = (tau/2) * G(c),   G(c) = sum_j (|y_j| - 1).
"""
import numpy as np


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley1(q):
    """Paley type I Hadamard matrix of order q+1, q = 3 mod 4 prime."""
    assert q % 4 == 3
    Q = np.array([[legendre(j - i, q) for j in range(q)] for i in range(q)], dtype=np.int64)
    n = q + 1
    S = np.zeros((n, n), dtype=np.int64)
    S[0, 1:] = 1
    S[1:, 0] = -1
    S[1:, 1:] = Q
    return np.eye(n, dtype=np.int64) + S


def sylvester(m):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < m:
        H = np.block([[H, H], [H, -H]])
    return H


def normalize(H):
    H = H * H[:, [0]]
    M = H.shape[0]
    assert (H.T @ H == M * np.eye(M, dtype=np.int64)).all()
    assert (H[:, 0] == 1).all()
    return H


def flipped_profile(H, b, assign):
    """b unflipped copies of every label a=1..M-1, plus one flipped copy per row x of label assign[x]
    (entry at row x negated).  Returns S (M x r) and the column labels and flip rows (-1 = unflipped)."""
    M = H.shape[0]
    cols, lab, frow = [], [], []
    for a in range(1, M):
        for _ in range(b):
            cols.append(H[:, a].copy()); lab.append(a); frow.append(-1)
    for x in range(M):
        a = assign[x]
        col = H[:, a].copy(); col[x] = -col[x]
        cols.append(col); lab.append(a); frow.append(x)
    S = np.array(cols, dtype=np.int64).T
    return S, np.array(lab), np.array(frow)


def Gval(S, c):
    y = S.T @ c
    return int(np.abs(y).sum()) - S.shape[1]


def F_of(H, c):
    T = H.T @ c
    assert T[0] == 0
    return int(np.abs(T[1:]).sum())
