"""Referee helpers for lower.md (independent of round5/checks/*).

Paley-I matrix, written from the text of lower.md Sec. 2 (not from paley.py):
  rows: x in F_q (index x) and infinity (index q);
  cols: constant column (index 0) and labels a in F_q (index a+1);
  H(x,a) = -chi(a-x) - [a=x],  H(inf,a) = 1,  constant column = 1.
All arithmetic is in Python ints / numpy int64 with tiny entries (exact)."""
import numpy as np


def isprime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return False
        f += 2
    return True


def chi(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley_H(q):
    assert isprime(q) and q % 4 == 3
    M = q + 1
    H = np.zeros((M, M), dtype=np.int64)
    H[:, 0] = 1
    for a in range(q):
        H[q, a + 1] = 1
        for x in range(q):
            H[x, a + 1] = -chi(a - x, q) - (1 if a == x else 0)
    G = H.T @ H
    assert (G == M * np.eye(M, dtype=np.int64)).all()
    return H


def T_of(H, c):
    T = c @ H
    assert T[0] == int(np.sum(c))
    return T[1:]


def Fval(H, c):
    return int(np.abs(T_of(H, c)).sum())


def in_E(H, c):
    """pair, +-column, or (s H_a + s' H_b)/2 (checked by definition, by brute force)."""
    M = H.shape[0]
    c = np.asarray(c)
    if int(np.abs(c).sum()) == 2 and int(np.abs(c).max()) == 1:
        return True
    for a in range(1, M):
        if (c == H[:, a]).all() or (c == -H[:, a]).all():
            return True
    # half-sums: support exactly M/2, entries +-1
    if int((c != 0).sum()) != M // 2 or int(np.abs(c).max()) != 1:
        return False
    T = T_of(H, c)
    nz = np.nonzero(T)[0]
    if len(nz) != 2:
        return False
    a, b = nz[0] + 1, nz[1] + 1
    for s in (1, -1):
        for t in (1, -1):
            g = s * H[:, a] + t * H[:, b]
            if (g % 2 == 0).all() and (c == g // 2).all():
                return True
    return False
