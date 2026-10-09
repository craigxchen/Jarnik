"""Independent builder (round6) of the one-flip, capacity-two Paley profile of round5/lower.md §3.1.

Paley-I matrix of order M = q0 + 1 (q0 = 3 mod 4 prime): rows x in F_q0 (index x) and infinity
(index q0); nonconstant labels a in F_q0;  H(x, a) = -chi(a - x) - [a == x],  H(inf, a) = 1.
Physical columns (a, k), k = 0..b-1, all equal to H_a before flipping.  Canonical assignment:
alpha(x) = x for finite x (flip column (x, 0), sibling (x, 1)); alpha(inf) = astar with flip column
(astar, 2) and sibling (astar, 3).  Rows a_x = (1 + S_x)/2 in {0,1}^r, r = b (M - 1)."""
import numpy as np


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def paley_H(q0):
    M = q0 + 1
    H = np.ones((M, q0), dtype=np.int64)          # nonconstant labels only
    for x in range(q0):
        for a in range(q0):
            H[x, a] = -legendre(a - x, q0) - (1 if a == x else 0)
    full = np.hstack([np.ones((M, 1), dtype=np.int64), H])
    assert (full.T @ full == M * np.eye(M, dtype=np.int64)).all()
    return H


def build(q0, b, astar=0):
    H = paley_H(q0)
    M = q0 + 1
    cols = [(a, k) for a in range(q0) for k in range(b)]
    idx = {c: j for j, c in enumerate(cols)}
    S = np.zeros((M, len(cols)), dtype=np.int64)
    for j, (a, k) in enumerate(cols):
        S[:, j] = H[:, a]
    flip, sib, alpha = {}, {}, {}
    for x in range(q0):
        flip[x], sib[x], alpha[x] = idx[(x, 0)], idx[(x, 1)], x
    flip[q0], sib[q0], alpha[q0] = idx[(astar, 2)], idx[(astar, 3)], astar
    for x in range(M):
        S[x, flip[x]] *= -1
    A = (1 + S) // 2
    return H, S, A, flip, sib, alpha
