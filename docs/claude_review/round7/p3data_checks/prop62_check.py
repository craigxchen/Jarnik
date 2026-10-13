"""Checks for Lemma 6.1 / Proposition 6.2 of p3data.md (separation demand of equal spacing).

E_xy = 2 log((M-1)/|i-j|) for points at positions i/(M-1), i = 0..M-1 (any order: E is a relabelling).
(1) max row sum of E  <=  2(1 + log 2)(M-1) - 4   (proved; checked here in floating point, M = 3..400)
(2) Sylvester Hadamard of order n = 2^ceil(log2 M), rows restricted to M rows, constant columns dropped:
    every pair of rows is separated by exactly n/2 of the remaining columns (exact integers, M = 3..130)
(3) W_cut(E) by the full LP for M = 6, 8, 10, 12 against 2 lambda_max(E), pi lambda_max(E)
    and the bound 2 pi (1 + log 2)(M-1) of Prop. 6.2(a).
(4) Lemma 6.1 additivity of margins and character margins under concatenation (random profiles, exact
    integer exponents, float weights; identity check).
"""
import math, itertools
import numpy as np
from qsolve import lp_full

ok_all = True
# (1)
worst = 0.0
for M in range(3, 401):
    n = M - 1
    rs = max(sum(2 * math.log(n / abs(i - j)) for j in range(M) if j != i) for i in range(M))
    bound = 2 * (1 + math.log(2)) * n - 4
    worst = max(worst, rs / bound)
    if rs > bound + 1e-9:
        ok_all = False; print('row-sum bound FAILS at M', M, rs, bound)
print('(1) max over M=3..400 of rowsum/bound = %.4f' % worst)

# (2)
def sylvester(k):
    H = np.array([[1]], dtype=np.int64)
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    return H
for M in range(3, 131):
    k = math.ceil(math.log2(M)); n = 2 ** k
    H = sylvester(k)[:M]
    keep = [j for j in range(n) if len(set(H[:, j])) > 1]
    S = H[:, keep]
    G = S @ S.T
    sep = (len(keep) - G) // 2
    off = sep[~np.eye(M, dtype=bool)]
    if not (off == n // 2).all():
        ok_all = False; print('Sylvester separation FAILS at M', M)
    if len(keep) > 2 * M - 2:
        ok_all = False; print('column count > 2M-2 at M', M)
print('(2) Sylvester restricted designs: every pair separated exactly n/2 times, r <= 2M-2, M = 3..130')

# (3)
for M in (6, 8, 10, 12):
    E = np.array([[0 if i == j else 2 * math.log((M - 1) / abs(i - j)) for j in range(M)] for i in range(M)])
    lam = float(np.linalg.eigvalsh(E)[-1])
    Wc = lp_full(E)[0]
    print('(3) M=%d: 2 lam_max(E) = %.3f  W_cut(E) = %.3f  pi lam_max(E) = %.3f  2pi(1+log2)(M-1) = %.3f' % (
        M, 2 * lam, Wc, math.pi * lam, 2 * math.pi * (1 + math.log(2)) * (M - 1)))
    if not (2 * lam - 1e-6 <= Wc <= math.pi * lam + 1e-6):
        ok_all = False

# (4)
rng = np.random.default_rng(5)
for trial in range(200):
    M = int(rng.integers(3, 9)); r1 = int(rng.integers(1, 7)); r2 = int(rng.integers(1, 7))
    A1 = rng.integers(0, 2, size=(M, r1)); A2 = rng.integers(0, 2, size=(M, r2))
    w1 = rng.uniform(1, 5, r1); w2 = rng.uniform(1, 5, r2)
    A = np.hstack([A1, A2]); w = np.concatenate([w1, w2])
    def mu(A, w):
        d = np.array([[np.abs(A[x] - A[y]) @ w for y in range(M)] for x in range(M)])
        return d - w.sum() / 2
    def sc(A, w, c):
        return np.abs(c @ A) @ w - w.sum() / 2
    m = mu(A, w) - mu(A1, w1) - mu(A2, w2)
    if np.abs(m[~np.eye(M, dtype=bool)]).max() > 1e-9:
        ok_all = False
    for _ in range(20):
        c = rng.integers(-3, 4, size=M); c[-1] -= c.sum()
        if abs(sc(A, w, c) - sc(A1, w1, c) - sc(A2, w2, c)) > 1e-9:
            ok_all = False
print('(4) additivity of mu and s_c under concatenation: 200 random profiles x 20 characters')
print('ALL PROP 6.2 CHECKS PASSED' if ok_all else 'SOME CHECK FAILED')
