"""Model boxes (independent re-implementation of the record's sdc/pal/syl pattern sets; conventions of round12/fake.md §1)."""
import itertools
import numpy as np

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def paley_I(q):
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

def syl_patterns(k):
    M = 2 ** k
    pts = list(itertools.product([0, 1], repeat=k))
    pats = set()
    for u in itertools.product([0, 1], repeat=k):
        if not any(u):
            continue
        S = [i for i, x in enumerate(pts) if sum(a * b for a, b in zip(u, x)) % 2 == 1]
        pats.add(canon(M, S))
    return M, sorted(pats)

def box(M, pats, L):
    d = M - 1
    # pattern P (omits row 0) -> indicator over coordinates 1..M-1 (index x-1)
    inds = [[1 if (x + 1) in P else 0 for x in range(d)] for P in pats]
    A = []
    for k in itertools.product([-1, 0, 1], repeat=d):
        ok = True
        for ind in inds:
            s = 0
            for a, b in zip(ind, k):
                if a:
                    s += b
            if abs(s) > L:
                ok = False
                break
        if ok:
            A.append(k)
    return A, inds
