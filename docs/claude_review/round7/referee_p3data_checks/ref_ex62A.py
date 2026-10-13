"""Referee: Example 6.2A. Build the Kronecker profile, check margins/W bound, and check whether the
Remark 4.2 tie occurs: sum over all cross pairs of two classes of mu_xy(w) == 0 for every w.
Then LP: max over weights w_j >= w0 of min over D=0 pairs mu_xy -> must be 0 (so no distinct-prime realisation)."""
import numpy as np, itertools
from scipy.optimize import linprog

def sylvester(n):
    H = np.array([[1]])
    while H.shape[0] < n: H = np.block([[H, H], [H, -H]])
    return H

def hadamard_normalised(k1):
    if k1 == 2: return np.array([[1, 1], [1, -1]])
    H = sylvester(k1)  # powers of 2 suffice for a test
    return H

def build(m, k1):
    nA = 1
    while nA < m: nA *= 2
    h = sylvester(nA)[:m]
    if k1 == 2:
        sig = np.array([[1], [-1]])
    else:
        sig = hadamard_normalised(k1)[:, 1:]
    rows = []; cls = []
    for c in range(m):
        for i in range(k1):
            rows.append(np.kron(h[c], sig[i])); cls.append(c)
    return np.array(rows), np.array(cls), nA

rng = np.random.default_rng(0)
for (m, k1) in [(2, 2), (3, 2), (5, 2), (3, 4), (5, 4), (2, 8), (6, 8)]:
    S, cls, nA = build(m, k1); M, r = S.shape; k = k1 - 1
    assert all(abs(S[:, j].sum()) < M for j in range(r))
    G = S @ S.T
    same = (cls[:, None] == cls[None, :]) & ~np.eye(M, dtype=bool)
    cross = cls[:, None] != cls[None, :]
    assert (G[same] == -nA).all() and (G[cross] == 0).all()
    # tie: sum over x in class 0, y in class 1 of s_x s_y  (vector over columns) == 0 ?
    v = sum(S[x]*S[y] for x in np.where(cls == 0)[0] for y in np.where(cls == 1)[0])
    tie = not v.any()
    # LP: maximise t s.t. mu_xy(w) >= t on cross pairs, w_j >= 1 (scale-free), sum w = r (normalisation)
    iu = [(x, y) for x in range(M) for y in range(x+1, M) if cross[x, y]]
    A = np.array([-0.5*S[x]*S[y] for x, y in iu])  # mu = A w
    # vars (w, t): max t -> min -t ; -A w + t <= 0 ; w >= 1
    c = np.zeros(r+1); c[-1] = -1
    Aub = np.hstack([-A, np.ones((len(iu), 1))])
    res = linprog(c, A_ub=Aub, b_ub=np.zeros(len(iu)), A_eq=np.hstack([np.ones((1, r)), np.zeros((1, 1))]), b_eq=[2*r],
                  bounds=[(1, None)]*r + [(None, None)], method='highs')
    # also complementary rows?
    comp = any((S[x] == -S[y]).all() for x in range(M) for y in range(x+1, M))
    print(f'm={m} k+1={k1} M={M} r={r}: Gram ok; identical cross-class cancellation={tie}; complementary rows={comp}; '
          f'max_w min_cross mu = {-res.fun:.3e}')
