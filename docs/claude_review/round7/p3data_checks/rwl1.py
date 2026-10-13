"""Reweighted-L1 heuristic for the prime-floor (semi-continuous) cut problem.

  minimise sum_s lam_s  s.t.  mu(lam) >= B,  lam_s in {0} u [wmin, inf).
Step 1: iteratively solve LPs  min sum_s c_s lam_s,  c_s = 1 + beta/(lam_s^prev + eps), which pushes
        weight onto few cuts (log-penalty surrogate of the fixed charge wmin per cut).
Step 2: on the resulting support T (plus the cuts of best reduced cost), solve the semi-continuous
        MILP exactly (sc_milp, restricted list) -- an upper bound for the full problem.
Only upper bounds come out of this file; every final profile is re-verified elsewhere.
"""
import numpy as np
from scipy.optimize import linprog
from qsolve import sc_milp
from qcore import pair_coeffs, best_cuts
from qcore import pair_index


def lp_weighted(B, masks, cost):
    M = B.shape[0]
    iu = pair_index(M)
    P, S = pair_coeffs(masks, M)
    res = linprog(cost, A_ub=-P, b_ub=-B[iu], bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    return res.x


def rwl1(B, masks, wmin, rounds=12, beta=None, eps=None):
    beta = wmin if beta is None else beta
    eps = 0.05 * wmin if eps is None else eps
    n = len(masks)
    lam = lp_weighted(B, masks, np.ones(n))
    hist = []
    for _ in range(rounds):
        cost = 1 + beta / (lam + eps)
        lam = lp_weighted(B, masks, cost)
        hist.append(int((lam > 1e-7).sum()))
    return lam, hist


def rwl1_sc(B, masks, wmin, nu=None, extra=150, time_limit=120, rounds=12, betas=(0.5, 1.0, 2.0)):
    """Run rwl1 for several beta multipliers; for each, restricted sc MILP on support + extra cuts."""
    best = None
    good = []
    if nu is not None:
        good, _ = best_cuts(nu, extra)
        good = list(good)
    for bm in betas:
        lam, hist = rwl1(B, masks, wmin, rounds=rounds, beta=bm * wmin)
        T = [masks[i] for i in range(len(masks)) if lam[i] > 1e-7]
        cand = sorted(set(T) | set(good))
        U = max(4 * float(lam.max()), 8 * wmin)
        r = sc_milp(B, cand, wmin, U=U, time_limit=time_limit)
        r['beta_mult'] = bm; r['rw_support'] = len(T)
        if r['W'] is not None and (best is None or r['W'] < best['W']):
            best = r
    return best


def nearest_menu(lam_s, logs, K):
    """indices of the K primes whose log is nearest to lam_s."""
    import bisect
    k = bisect.bisect_left(logs, lam_s)
    lo, hi, out = k - 1, k, []
    while len(out) < K and (lo >= 0 or hi < len(logs)):
        if hi >= len(logs) or (lo >= 0 and lam_s - logs[lo] <= logs[hi] - lam_s):
            out.append(lo); lo -= 1
        else:
            out.append(hi); hi += 1
    return out


def menu_milp(B, masks, lam, plist, K=10, time_limit=120):
    """every support cut takes exactly one prime among the K nearest (in log) to its weight; primes
    distinct; all pair margins >= B + 1e-7; minimise total log weight.  Returns cols or None."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import csr_matrix
    import math
    M = B.shape[0]
    iu = pair_index(M)
    logs = [math.log(p) for p in plist]
    T = [(m, l) for m, l in zip(masks, lam) if l > 0]
    cuts = [m for m, _ in T]
    opts = [nearest_menu(l, logs, K) for _, l in T]
    P, S = pair_coeffs(cuts, M)
    var = [(c, k) for c in range(len(cuts)) for k in opts[c]]
    nv = len(var)
    cost = np.array([logs[k] for c, k in var])
    Ap = np.zeros((len(iu[0]), nv)); Ac = np.zeros((len(cuts), nv))
    used = sorted({k for _, k in var}); pos = {k: i for i, k in enumerate(used)}
    Au = np.zeros((len(used), nv))
    for v, (c, k) in enumerate(var):
        Ap[:, v] = P[:, c] * logs[k]; Ac[c, v] = 1; Au[pos[k], v] = 1
    A = csr_matrix(np.vstack([Ap, Ac, Au]))
    lo = np.concatenate([B[iu] + 1e-7, np.ones(len(cuts)), np.full(len(used), -np.inf)])
    hi = np.concatenate([np.full(len(iu[0]), np.inf), np.ones(len(cuts)), np.ones(len(used))])
    res = milp(cost, constraints=[LinearConstraint(A, lo, hi)], integrality=np.ones(nv),
               bounds=Bounds(np.zeros(nv), np.ones(nv)), options=dict(time_limit=time_limit, disp=False))
    if res.x is None:
        return None
    cols = []
    for v, (c, k) in enumerate(var):
        if res.x[v] > 0.5:
            cols.append((cuts[c], [plist[k]]))
    return cols
