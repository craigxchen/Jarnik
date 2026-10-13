"""Flip local search (FLS) for the prime-floor cut problem with a fixed number r of columns.

Columns are +-1 vectors S[:, j] (nonconstant), weights w_j >= w0.  For fixed S the best weights solve
    LP(S):  min sum_j w_j  s.t.  sum_j w_j ([S_xj != S_yj] - 1/2) >= B_xy (x<y),  w_j >= w0,
with an l1 penalty on violations (phase-1 style) so that every S has a value.  Local moves: flip one
entry S_xj (columns must stay nonconstant), ranked by the first-order gain from the LP duals and
accepted if the re-solved LP value drops.  Output: upper bounds only (every final solution is a
prime-floor-feasible weight vector; the realisation with actual primes is done elsewhere).
"""
import numpy as np
from scipy.optimize import linprog
from qcore import pair_index

PEN = 1e3


def lp_S(B, S, w0):
    M, r = S.shape
    iu = pair_index(M)
    sep = (S[iu[0], :] != S[iu[1], :]).astype(float)          # npair x r
    P = sep - 0.5
    npair = P.shape[0]
    # variables: w (r), slack (npair)
    c = np.concatenate([np.ones(r), PEN * np.ones(npair)])
    A = np.hstack([-P, -np.eye(npair)])
    b = -B[iu]
    bounds = [(w0, None)] * r + [(0, None)] * npair
    res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
    if res.status != 0:
        return np.inf, None, None, None
    w = res.x[:r]; sl = res.x[r:]
    nu = -res.ineqlin.marginals
    return res.fun, w, sl.sum(), nu


def fls(B, S0, w0, max_rounds=200, cand=12, rng=None, verbose=False):
    rng = np.random.default_rng(0) if rng is None else rng
    S = S0.copy()
    M, r = S.shape
    iu = pair_index(M)
    val, w, sl, nu = lp_S(B, S, w0)
    for rnd in range(max_rounds):
        N = np.zeros((M, M)); N[iu] = nu; N = N + N.T
        # gain of flipping (x, j): margins of pairs (x, y) change by w_j (1 - 2 sep_xyj)
        same = (S[:, None, :] == S[None, :, :]).astype(float)   # M x M x r
        delta = (2 * same - 1) * w[None, None, :]                # change of mu_xy if x flips in column j
        gain = np.einsum('xy,xyj->xj', N, delta)                  # first-order decrease of the LP value
        order = np.dstack(np.unravel_index(np.argsort(-gain, axis=None), gain.shape))[0]
        improved = False
        tried = 0
        for x, j in order:
            if gain[x, j] <= 1e-9 or tried >= cand:
                break
            col = S[:, j].copy(); col[x] = -col[x]
            if abs(col.sum()) == M:
                continue
            tried += 1
            S2 = S.copy(); S2[:, j] = col
            v2, w2, sl2, nu2 = lp_S(B, S2, w0)
            if v2 < val - 1e-7:
                S, val, w, sl, nu = S2, v2, w2, sl2, nu2
                improved = True
                break
        if not improved:
            # random perturbation of one entry among the top gains, accept only if not worse by > 2%
            break
        if verbose and rnd % 20 == 0:
            print('  fls round %d value %.3f slack %.3g' % (rnd, val, sl))
    return dict(W=val, w=w, S=S, slack=sl)


def sylvester_rows(M):
    k = int(np.ceil(np.log2(M)))
    H = np.array([[1]])
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    H = H[:M]
    keep = [j for j in range(H.shape[1]) if abs(H[:, j].sum()) < M]
    return H[:, keep]


def init_from_cuts(masks, lam, M, r):
    from qcore import mask_to_s
    idx = np.argsort(-np.asarray(lam))[:r]
    return np.array([mask_to_s(masks[i], M) for i in idx]).T


def best_fls(B, w0, inits, rounds=200, cand=12):
    best = None
    for name, S0 in inits:
        out = fls(B, S0, w0, max_rounds=rounds, cand=cand)
        out['init'] = name; out['r'] = S0.shape[1]
        if out['slack'] < 1e-7 and (best is None or out['W'] < best['W']):
            best = out
    return best


def drop_columns(B, S, w0):
    """greedily drop columns while the LP value (with zero slack) decreases."""
    val, w, sl, nu = lp_S(B, S, w0)
    improved = True
    while improved and S.shape[1] > 1:
        improved = False
        for j in np.argsort(w):
            S2 = np.delete(S, j, axis=1)
            v2, w2, sl2, _ = lp_S(B, S2, w0)
            if sl2 is not None and sl2 < 1e-7 and v2 < val - 1e-7:
                S, val, w, sl = S2, v2, w2, sl2
                improved = True
                break
    return S, val


def ils(B, S0, w0, kicks=60, kick_size=3, cand=40, rng=None):
    """iterated local search: fls, then repeated random kicks of the best solution."""
    rng = np.random.default_rng(0) if rng is None else rng
    best = fls(B, S0, w0, cand=cand, rng=rng)
    M = S0.shape[0]
    for _ in range(kicks):
        S = best['S'].copy()
        r = S.shape[1]
        for _ in range(kick_size):
            x, j = rng.integers(M), rng.integers(r)
            S[x, j] = -S[x, j]
            if abs(S[:, j].sum()) == M:
                S[x, j] = -S[x, j]
        out = fls(B, S, w0, cand=cand, rng=rng)
        if out['W'] < best['W'] - 1e-7:
            best = out
    if best['slack'] < 1e-7:
        S, val = drop_columns(B, best['S'], w0)
        if val < best['W']:
            v, w, sl, _ = lp_S(B, S, w0)
            best = dict(W=v, w=w, S=S, slack=sl)
    return best


def best_ils(B, w0, inits, kicks=60, cand=40, seed=0):
    rng = np.random.default_rng(seed)
    best = None
    for name, S0 in inits:
        out = ils(B, S0, w0, kicks=kicks, cand=cand, rng=rng)
        out['init'] = name; out['r'] = out['S'].shape[1]
        if out['slack'] < 1e-7 and (best is None or out['W'] < best['W']):
            best = out
    return best
