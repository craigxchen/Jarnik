"""[J] Half-row Paley pattern (meaning referee item 3): S = K - I, K_iv = chi(v - i), m = 3 mod 4 prime;
retained rows R = {i : S_i0 = +1}, missing rows M = the rest; T = [ones; S_M].
X.3's deletion argument excludes R (for P >= P0(m, C), two primes per column) iff ker T is spanned by
integer vectors of support < m/2 (strictly: sum of log Q_v over the support below log sqrt q - O(log mC)).
Compute the minimum support (spark) of T and whether vectors of support <= floor((m-1)/2) span ker T."""
import sys
from itertools import combinations
import numpy as np
from sympy import Matrix


def chi(a, m):
    a %= m
    if a == 0:
        return 0
    return 1 if pow(a, (m - 1) // 2, m) == 1 else -1


for m in [int(x) for x in sys.argv[1:]] or [11, 19]:
    S = np.array([[(-1 if i == v else chi(v - i, m)) for v in range(m)] for i in range(m)])
    assert (S @ S.T == (m + 1) * np.eye(m, dtype=int) - 1).all()
    R = [i for i in range(m) if S[i, 0] == 1]
    Mrows = [i for i in range(m) if S[i, 0] != 1]
    T = np.vstack([np.ones(m, dtype=int), S[Mrows]])
    rT = np.linalg.matrix_rank(T)
    dimker = m - rT
    smax = (m - 1) // 2
    found = {}
    for k in range(2, smax + 1):
        vecs = []
        for J in combinations(range(m), k):
            sub = T[:, J]
            if np.linalg.matrix_rank(sub) < k:
                ns = Matrix(sub.tolist()).nullspace()
                for v in ns:
                    full = [0] * m
                    for j, x in zip(J, v):
                        full[j] = x
                    vecs.append(full)
        if vecs:
            found[k] = vecs
        allv = [v for kk in found for v in found[kk]]
        span = Matrix(allv).rank() if allv else 0
        print("m=%d retained=%d T=%dx%d dim ker=%d: support<=%d vectors found=%d, span=%d"
              % (m, len(R), T.shape[0], m, dimker, k, len(allv), span), flush=True)
        if span == dimker:
            break
