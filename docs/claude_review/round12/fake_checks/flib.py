"""Library for the round-12 fake test (T = 1 slab-boxes of barrier profiles).

Conventions.  Rows (points) 0..M-1.  A sign column is recorded as a pattern S (set of rows with entry -1 after
normalisation); S and its complement give the same slab.  Base point b = 0 (row 0); coordinates of k' are rows
1..M-1 (index x-1).  For a pattern S, relative pattern P = S if b not in S else complement(S); then
<k, beta> = +-2 k(P) on Z^M_0 and the slab |<k,beta>| <= 2L reads |k(P)| <= L.
T = 1 box: A = {k' in {-1,0,1}^d : |k(P)| <= L for all patterns}.
"""
import itertools
import numpy as np


def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)


def paley_I(q):
    """Paley-I Hadamard matrix of order q+1 (q = 3 mod 4), rows/cols indexed [inf] + F_q (as in round11 model8.py)."""
    idx = ['inf'] + list(range(q))
    M = q + 1
    H = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf':
                H[a, b] = 1
            elif y == 'inf':
                H[a, b] = -1
            else:
                H[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    return H


def canon(M, S):
    S = frozenset(S)
    if 0 in S:
        S = frozenset(range(M)) - S
    return tuple(sorted(S))


def sdc_patterns(q):
    """Skew double core (round5/exact.md Thm 4.1) at M = q+1: cuts S_a, S'_a = S_a xor {0,a}."""
    K = paley_I(q)
    M = q + 1
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all():
        K = K.T
    assert ((K + K.T) == 2 * np.eye(M, dtype=int)).all()
    Dg = np.diag(K[:, 0])
    K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        Sp = sorted(set(S) ^ {0, a})
        pats.add(canon(M, S))
        pats.add(canon(M, Sp))
    return M, sorted(pats)


def pal_patterns(q):
    """Round-5 Paley-I fake profile P(H, pi, b, a*) (construct.md §2), distinct patterns only (copies merge at T=1):
    unflipped labels a = 1..M-1 and for every row x the label a(x) flipped in row x (pi(i) = i, a(inf) = a* = 1)."""
    H = paley_I(q)
    M = q + 1
    Hn = H.copy()
    Hn[1:, :] *= -1
    lab = {0: 1}
    for i in range(q):
        lab[1 + i] = 1 + i
    pats = set()
    for a in range(1, M):
        S = tuple(x for x in range(M) if Hn[x, a] == 1)
        if 0 < len(S) < M:
            pats.add(canon(M, S))
    for x in range(M):
        c = Hn[:, lab[x]].copy()
        c[x] *= -1
        S = tuple(y for y in range(M) if c[y] == 1)
        if 0 < len(S) < M:
            pats.add(canon(M, S))
    return M, sorted(pats)


def sylvester_patterns(k):
    M = 2 ** k
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if bin(a & x).count('1') % 2 == 1]
        pats.add(canon(M, S))
    return M, sorted(pats)


def rel(P):
    """pattern (rows, base 0 excluded) -> coordinate indices (row x -> x-1)."""
    return [x - 1 for x in P]


def box(M, pats, L):
    """Enumerate A by depth-first search with pruning on partial slab sums."""
    d = M - 1
    cols = [rel(P) for P in pats]
    # for each coordinate, list of columns containing it; remaining counts for pruning
    inc = [[j for j, C in enumerate(cols) if x in C] for x in range(d)]
    rem = [[sum(1 for y in C if y >= x) for C in cols] for x in range(d + 1)]
    out = []
    k = [0] * d
    sums = [0] * len(cols)

    def rec(x):
        if x == d:
            out.append(tuple(k))
            return
        for v in (-1, 0, 1):
            ok = True
            for j in inc[x]:
                sums[j] += v
            # prune: |sum_j| - remaining(after x) > L
            for j in range(len(cols)):
                if abs(sums[j]) - rem[x + 1][j] > L:
                    ok = False
                    break
            if ok:
                k[x] = v
                rec(x + 1)
            for j in inc[x]:
                sums[j] -= v
        k[x] = 0

    rec(0)
    return out


def reduced_monomials(d, deg):
    """exponent vectors a in {0,1,2}^d with |a| = deg."""
    res = []

    def rec(x, left, cur):
        if x == d:
            if left == 0:
                res.append(tuple(cur))
            return
        if left > 2 * (d - x):
            return
        for e in (0, 1, 2):
            if e <= left:
                cur.append(e)
                rec(x + 1, left - e, cur)
                cur.pop()

    rec(0, deg, [])
    return res
