"""Referee's own library: sdc_q from round5/exact.md Thm 4.1 with the Paley-I skew-Hadamard matrix.
Index 0 = inf, index 1+x = x in F_q."""
import numpy as np

def chi_table(q):
    t = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        t[x] = 1 if pow(x, (q-1)//2, q) == 1 else -1
    return t

def skew_hadamard(q):
    ch = chi_table(q); M = q+1
    K = np.zeros((M, M), dtype=np.int64)
    K[0, :] = 1; K[:, 0] = -1; K[0, 0] = 1
    for x in range(q):
        for y in range(q):
            K[1+x, 1+y] = (1 if x == y else 0) + ch[(y-x) % q]
    assert ((K + K.T) == 2*np.eye(M, dtype=np.int64)).all()
    assert (K @ K.T == M*np.eye(M, dtype=np.int64)).all()
    D = np.diag(K[:, 0]); K = D @ K @ D
    assert (K[:, 0] == 1).all()
    return K

def patterns(q):
    """incidence matrix (2q x (q+1)) of cuts S_a, S'_a = S_a xor {0,a} (exact.md Thm 4.1)."""
    K = skew_hadamard(q); M = q+1
    rows = []
    for a in range(1, M):
        S = (K[:, a] == -1).astype(np.int64)
        Sp = S.copy(); Sp[0] ^= 1; Sp[a] ^= 1
        rows.append(S); rows.append(Sp)
    P = np.array(rows)
    assert len({tuple(r) for r in P} | {tuple(1-r) for r in P}) == 4*q  # distinct as partitions
    assert (P.sum(1) == M//2).all()
    return P

def sigma_def(P, lam):
    """sigma = sum_S (|lam(S)| - 1/2), exact integers -> returns 2*sigma as int to stay integral."""
    v = P @ lam
    return int(2*np.abs(v).sum() - P.shape[0])  # = 2 sigma
