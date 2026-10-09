"""Shared helpers: the normalised Paley-I matrix of order M = q+1 (q = 3 mod 4 prime), in the
convention of round4/flipped.md Lemma 6.2:

  rows    x in F_q  (indices 0..q-1)  and  infinity (index q)
  columns the constant column (index 0)  and  a in F_q (index a+1)
  H(x, a) = -chi(a - x) - [a == x],   H(inf, a) = 1,   constant column = 1.

All arithmetic is exact (numpy int64 on small integers, or Python ints)."""
import numpy as np


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley(q):
    assert is_prime(q) and q % 4 == 3
    M = q + 1
    H = np.ones((M, M), dtype=np.int64)
    for x in range(q):
        for a in range(q):
            H[x, a + 1] = -legendre(a - x, q) - (1 if a == x else 0)
    # row q = infinity: all ones (already); column 0 constant ones (already)
    assert (H.T @ H == M * np.eye(M, dtype=np.int64)).all(), "not Hadamard"
    assert (H[:, 0] == 1).all()
    return H


def transform(H, c):
    """T = H_A^T c over the nonconstant labels (the constant label gives sum c = 0)."""
    T = c @ H
    assert T[0] == 0
    return T[1:]


def F(H, c):
    return int(np.abs(transform(H, c)).sum())


def exceptional(H, c):
    """pair, +-column, or +-(H_a +- H_b)/2 ?"""
    M = H.shape[0]
    n = int(np.abs(c).sum())
    if n == 2:
        return "pair"
    T = transform(H, c)
    nz = np.nonzero(T)[0]
    if len(nz) == 1 and abs(T[nz[0]]) == M:
        a = nz[0] + 1
        if (c == H[:, a]).all() or (c == -H[:, a]).all():
            return "column"
    if len(nz) == 2 and abs(T[nz[0]]) == M // 2 and abs(T[nz[1]]) == M // 2:
        a, b = nz[0] + 1, nz[1] + 1
        for s1 in (1, -1):
            for s2 in (1, -1):
                g = (s1 * H[:, a] + s2 * H[:, b])
                if (g % 2 == 0).all() and (c == g // 2).all():
                    return "halfsum"
    return None
